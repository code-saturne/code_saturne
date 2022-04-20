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
#include "bft_printf.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the solution has a correct size and if needed allocate and
 *         initialize a new solution array before the solving step
 *
 *         - The scatter view corresponds to the (partitioned) mesh view
 *         - The gather view is purely algebraic => matrix view.
 *         - In parallel computation, n_cols > n_rows (there are elements to
 *           consider through the parallel interfaces)
 *
 * \param[in]  stride         number of DoFs by element
 * \param[in]  n_scatter_elts local number of elements (may be != n_gather_elts)
 * \param[in]  x              current pointer to the solution array
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 *
 * \return the pointer to the solution array to use
 */
/*----------------------------------------------------------------------------*/

static cs_real_t  *
_set_xsol(int                 stride,
          cs_lnum_t           n_scatter_elts,
          cs_real_t          *x,
          const cs_matrix_t  *matrix)
{
  cs_real_t  *xsol = NULL;

  /* The number of rows in the matrix is equal to the number of elements in the
     "gather" view: n_rows = n_gather_elts */

  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);

#if defined(DEBUG) && !defined(NDEBUG)
  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(matrix);
  assert(n_cols >= n_rows);
  assert(n_rows <= n_scatter_elts);
#endif

  if (n_cols > n_scatter_elts) {

    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, stride*n_cols, cs_real_t);
    memcpy(xsol, x, stride*n_scatter_elts*sizeof(cs_real_t));

  }
  else
    xsol = x;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 0
  bft_printf(" %s >> rank_id:%d\n"
             " stride:           %d\n"
             " n_scatter_elts:   %d\n"
             " n_matrix_rows:    %d\n"
             " n_matrix_columns: %d\n",
             __func__, cs_glob_rank_id, stride, n_scatter_elts, n_rows, n_cols);
#endif

  return xsol;
}

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
 * \param[in]      interlace  is data interlaced or not
 * \param[in]      x_size     size of the vector unknowns (scatter view)
 * \param[in]      rset       pointer to a range set structure
 * \param[in]      rhs_redux  do or not a parallel sum reduction on the RHS
 * \param[in, out] x          array of unknowns (in: initial guess)
 * \param[in, out] b          right-hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_solve_prepare_system(int                     stride,
                            bool                    interlace,
                            cs_lnum_t               x_size,
                            const cs_range_set_t   *rset,
                            bool                    rhs_redux,
                            cs_real_t              *x,
                            cs_real_t              *b)
{
  if (rset == NULL)
    return;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 0
  bft_printf(" %s >> rank_id:%d\n"
             " stride:           %d\n"
             " interlaced:       %s\n"
             " rset->n_elts      %d %d %d\n"
             " rset->l_range     %u %u\n",
             __func__, cs_glob_rank_id, stride, cs_base_strtf(interlace),
             rset->n_elts[0], rset->n_elts[1], rset->n_elts[2],
             rset->l_range[0], rset->l_range[1]);
#endif

  /* x and b should be changed to have a "gathered" view through the range set
     operation.  Their size is equal to n_sles_gather_elts <=
     n_sles_scatter_elts */

  /* Compact numbering to fit the algebraic decomposition */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 2
  cs_dbg_darray_to_listing(" Initial guess (scatter view)",
                           x_size*stride, x, 4*stride);
#endif

  cs_range_set_gather(rset,
                      CS_REAL_TYPE, /* type */
                      stride,       /* stride */
                      x,            /* in: size = n_sles_scatter_elts */
                      x);           /* out: size = n_sles_gather_elts */

  /* The right-hand side stems from a cellwise building on this rank.
     Other contributions from distant ranks may contribute to an element
     owned by the local rank */

  if (rhs_redux && rset->ifs != NULL)
    cs_interface_set_sum(rset->ifs,
                         x_size, stride, interlace, CS_REAL_TYPE,
                         b);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 2
  cs_dbg_darray_to_listing(" RHS (scatter view)",
                           x_size*stride, b, 4*stride);
#endif

  cs_range_set_gather(rset,
                      CS_REAL_TYPE,/* type */
                      stride,      /* stride */
                      b,           /* in: size = n_sles_scatter_elts */
                      b);          /* out: size = n_sles_gather_elts */
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
  /* Set xsol (manage allocation and initialization in case of parallelism) */

  cs_real_t  *xsol = _set_xsol(1, n_scatter_dofs, x, matrix);

  /* Prepare solving (handle parallelism) and switch to a "gather" view
   * stride = 1 for scalar-valued */

  cs_cdo_solve_prepare_system(1, false, n_scatter_dofs, rset, rhs_redux,
                              xsol, b);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 2
  if (slesp->verbosity > 2)
    cs_dbg_dump_local_scalar_msr_matrix(slesp->name, matrix);
#endif

  /* Retrieve the solving info structure stored in the cs_field_t structure
     and then solve the linear system */

  cs_field_t  *fld = cs_field_by_id(slesp->field_id);
  cs_solving_info_t  sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = normalization;

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

  /* Switch back from the "gather" view to the "scatter" view */

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

  if (xsol != x)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with vector-valued
 *         degrees of freedom (DoFs).
 *         Number of DoFs is equal to 3*n_scatter_elts
 *
 * \param[in]  n_scatter_elts local number of elements (may be != n_gather_elts)
 * \param[in]  interlace      way to arrange data (true/false)
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  rset           pointer to a cs_range_set_t structure
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
cs_cdo_solve_vector_system(cs_lnum_t                     n_scatter_elts,
                           bool                          interlace,
                           const cs_param_sles_t        *slesp,
                           const cs_matrix_t            *matrix,
                           const cs_range_set_t         *rset,
                           cs_real_t                     normalization,
                           bool                          rhs_redux,
                           cs_sles_t                    *sles,
                           cs_real_t                    *x,
                           cs_real_t                    *b)
{
  /* Set xsol (manage allocation and initialization in case of parallelism) */

  cs_real_t  *xsol = _set_xsol(3, n_scatter_elts, x, matrix);

  /* Prepare solving (handle parallelism) and switch to a "gather" view */

  cs_cdo_solve_prepare_system(3, interlace, n_scatter_elts, rset, rhs_redux,
                              xsol, b);

  /* Retrieve the solving info structure stored in the cs_field_t structure
     and then solve the linear system */

  cs_field_t  *fld = cs_field_by_id(slesp->field_id);
  cs_solving_info_t  sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = normalization;

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

  /* Switch back from the "gather" view to the "scatter" view */

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 3, /* type and stride */
                       xsol, x);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 3, /* type and stride */
                       b, b);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 1
  cs_dbg_fprintf_system(slesp->name, cs_cdo_solve_dbg_counter++,
                        slesp->verbosity,
                        x, b, 3*n_scatter_elts);
#endif

  if (xsol != x)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
