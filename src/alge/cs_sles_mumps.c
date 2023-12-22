/*============================================================================
 * Sparse Linear Equation Solvers using MUMPS (a sparse direct solver library)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <sys/sysinfo.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * MUMPS headers
 *----------------------------------------------------------------------------*/

#include <dmumps_c.h>
#include <smumps_c.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_base.h"
#include "cs_log.h"
#include "cs_fp_exception.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_parall.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_sles_mumps.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_mumps.c

  \brief Set of functions to handle the interface with the MUMPS library
         to solve sparse linear system with a direct approach

  \page sles_mumps MUMPS-based linear solvers.

  \typedef cs_sles_mumps_setup_hook_t

  \brief Function pointer for user settings of a MUMPS DMUMPS_STRUC_C
  solver setup.

  This function is called the end of the setup stage.

  Note: if the context pointer is non-NULL, it must point to valid data
  when the selection function is called so that value or structure should
  not be temporary (i.e. local);

  \param[in, out]  context  pointer to optional (untyped) value or structure
  \param[in, out]  dmumps   pointer to DMUMPS_STRUC_C structure

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* MUMPS code to detect the calculation step */

#define MUMPS_JOB_END             -2
#define MUMPS_JOB_INIT            -1
#define MUMPS_JOB_ANALYSIS         1
#define MUMPS_JOB_FACTORIZATION    2
#define MUMPS_JOB_SOLVE            3

/* Default number for MPI_COMM_WORLD */

#define USE_COMM_WORLD         -987654

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

static const int  cs_sles_mumps_n_max_tries = 3;

/* Type of factorization to perform with MUMPS (precision / facto) */

typedef enum {

  CS_SLES_MUMPS_DOUBLE_LU,       /* LU facto. with dmumps */
  CS_SLES_MUMPS_DOUBLE_LDLT_SYM, /* LDLt facto. with dmumps for sym. matrices */
  CS_SLES_MUMPS_DOUBLE_LDLT_SPD, /* LDLt facto. with dmumps for SPD matrices */
  CS_SLES_MUMPS_SINGLE_LU,       /* LU facto. with smumps */
  CS_SLES_MUMPS_SINGLE_LDLT_SYM, /* LDLt facto. with smumps for sym. matrices */
  CS_SLES_MUMPS_SINGLE_LDLT_SPD, /* LDLt facto. with smumps for SPD matrices */

  CS_SLES_MUMPS_N_TYPES

} cs_sles_mumps_type_t;


struct _cs_sles_mumps_t {

  cs_sles_mumps_type_t    type;        /* Type of usage of MUMPS */

  /* Performance data */

  int                  n_tries;       /* Number of analysis/facto done */
  int                  n_setups;      /* Number of times system setup */
  int                  n_solves;      /* Number of times system solved since
                                       * it's a direct solver thus
                                       * n_solves = n_iterations_tot */

  cs_timer_counter_t   t_setup;       /* Total setup (factorization) */
  cs_timer_counter_t   t_solve;       /* Total time used */

  /* Additional setup options when used */

  bool                         is_pc;         /* MUMPS is used as a
                                                 preconditioner */
  const cs_matrix_t           *matrix;        /* Only useful when used as
                                                 preconditioner */

  const cs_param_sles_t       *sles_param;    /* set of parameter for SLES */

  void                        *hook_context;  /* Optional user context */
  cs_sles_mumps_setup_hook_t  *setup_hook;    /* Post setup function */

  /* MUMPS structure: either a pointer to SMUMPS_STRUC_C structure in case of a
   * single-precision settings or a pointer to a DMUMPS_STRUC_C structure in
   * case of a double-precision settings */

  void                        *mumps_struct;

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static int  _n_mumps_systems = 0;
static double  cs_sles_mumps_zero_dthreshold = 64*DBL_MIN;
static double  cs_sles_mumps_zero_fthreshold = 128*FLT_MIN;

/*============================================================================
 * Static inline private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a periodicity has to be performed during the computation
 *
 * \param[in] halo   pointer to a halo structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_have_perio(const cs_halo_t     *halo)
{
  bool  have_perio = false;
  if (halo != NULL)
    if (halo->n_transforms > 0)
      have_perio = true;

  if (have_perio)
    bft_error(__FILE__, __LINE__, 0,
              " %s: No periodicity requested with MUMPS up to now.\n",
              __func__);

  return have_perio;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the MUMPS solver is defined as a double-precision or a
 *        single-precision solver and for which kind of factorization
 *
 * \param[in] slesp   pointer to a SLES parameter structure
 *
 * \return the type of MUMPS usage
 */
/*----------------------------------------------------------------------------*/

static inline cs_sles_mumps_type_t
_set_type(const cs_param_sles_t  *slesp)
{
  assert(slesp != NULL);
  assert(slesp->context_param != NULL);

  cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

  if (mumpsp->is_single) {

    switch(mumpsp->facto_type) {

    case CS_PARAM_SLES_FACTO_LU:
      return CS_SLES_MUMPS_SINGLE_LU;
    case CS_PARAM_SLES_FACTO_LDLT_SYM:
      return CS_SLES_MUMPS_SINGLE_LDLT_SYM;
    case CS_PARAM_SLES_FACTO_LDLT_SPD:
      return CS_SLES_MUMPS_SINGLE_LDLT_SPD;

    default:
      return CS_SLES_MUMPS_N_TYPES;
    }

  }
  else {

    switch(mumpsp->facto_type) {

    case CS_PARAM_SLES_FACTO_LU:
      return CS_SLES_MUMPS_DOUBLE_LU;
    case CS_PARAM_SLES_FACTO_LDLT_SYM:
      return CS_SLES_MUMPS_DOUBLE_LDLT_SYM;
    case CS_PARAM_SLES_FACTO_LDLT_SPD:
      return CS_SLES_MUMPS_DOUBLE_LDLT_SPD;

    default:
      return CS_SLES_MUMPS_N_TYPES;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if MUMPS is used as preconditioner
 *
 * \param[in] slesp   pointer to a SLES parameter structure
 *
 * \return true if MUMPS is used as preconditioner, false otherwise
 */
/*----------------------------------------------------------------------------*/

static inline bool
_set_pc_usage(const cs_param_sles_t  *slesp)
{
  assert(slesp != NULL);

  switch (slesp->solver) {

  case CS_PARAM_ITSOL_MUMPS:
    return false;

  default: /* Not a solver. Try as preconditioner */
    switch(slesp->precond) {

    case CS_PARAM_PRECOND_MUMPS:
      return true;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: MUMPS not defined as solver or as preconditioner\n"
                "%s: for the system \"%s\".\n",
                __func__, __func__, slesp->name);
      break;
    }

  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the MUMPS solver is defined as a double-precision or a
 *        single-precision solver.
 *
 * \param[in] c   pointer to the context structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_is_dmumps(const cs_sles_mumps_t  *c)
{
  switch (c->type) {

  case CS_SLES_MUMPS_DOUBLE_LDLT_SPD:
  case CS_SLES_MUMPS_DOUBLE_LDLT_SYM:
  case CS_SLES_MUMPS_DOUBLE_LU:
    return true;

  case CS_SLES_MUMPS_SINGLE_LDLT_SPD:
  case CS_SLES_MUMPS_SINGLE_LDLT_SYM:
  case CS_SLES_MUMPS_SINGLE_LU:
    return false;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Undefined MUMPS type for the system \"%s\".\n",
              __func__, c->sles_param->name);
    break;
  }

  return false;
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the MUMPS structure in case of double-precision
 *        computation. Default settings.
 *
 * \param[in] context   pointer to the preconditioner context
 * \param[in] logging   if true, logging description; if false, canonical name
 */
/*----------------------------------------------------------------------------*/

static const char *
_mumps_pc_get_type(const void  *context,
                   bool         logging)
{
  CS_UNUSED(context);

  if (logging == false) {
    static const char t[] = "mumps preconditioner";
    return t;
  }
  else {
    static const char t[] = N_("MUMPS preconditioner");
    return _(t);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for setup of MUMPS as preconditioner
 *
 * \param[in, out] context    pointer to preconditioner context
 * \param[in]      name       pointer to name of associated linear system
 * \param[in]      a          matrix
 * \param[in]      accel      use accelerator version ?
 * \param[in]      verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

static void
_mumps_pc_setup(void               *context,
                const char         *name,
                const cs_matrix_t  *a,
                bool                accel,
                int                 verbosity)
{
  CS_UNUSED(accel);

  cs_sles_mumps_t  *c = context;

  c->matrix = a;                /* Only a shared pointer */

  cs_sles_mumps_setup(context, name, a, verbosity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply MUMPS as preconditioner
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contains the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).

 * \param[in, out] context    pointer to preconditioner context
 * \param[in]      x_in       input_vector
 * \param[in,out]  x_out      input/output vector
 *
 * \return the preconditioner status
 */
/*----------------------------------------------------------------------------*/

static cs_sles_pc_state_t
_mumps_pc_apply(void                *context,
                const cs_real_t     *x_in,
                cs_real_t           *x_out)
{
  int     n_iter;
  double  residual;
  cs_real_t  *_rhs = NULL;

  cs_sles_mumps_t  *c = context;
  assert(c->is_pc && c->matrix != NULL);

  const cs_real_t *rhs = x_in;
  const cs_param_sles_t  *slesp = c->sles_param;
  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(c->matrix);
  const cs_lnum_t n_rows = cs_matrix_get_n_rows(c->matrix) * db_size;

  /* If preconditioner is "in-place", use additional buffer */

  if (x_in == NULL) {

    const cs_lnum_t n_cols = cs_matrix_get_n_columns(c->matrix) * db_size;
    BFT_MALLOC(_rhs, n_cols, cs_real_t);

    cs_array_real_copy(n_rows, x_out, _rhs);
    rhs = _rhs;

  }

  cs_array_real_fill_zero(n_rows, x_out);

  cs_sles_convergence_state_t  cvg = cs_sles_mumps_solve(context,
                                                         slesp->name,
                                                         c->matrix,
                                                         slesp->verbosity,
                                                         slesp->cvg_param.rtol,
                                                         1.0,
                                                         &n_iter,
                                                         &residual,
                                                         rhs,
                                                         x_out,
                                                         0,
                                                         NULL);

  if (x_in == NULL)
    BFT_FREE(_rhs);

  cs_sles_pc_state_t state;

  switch(cvg) {
  case CS_SLES_DIVERGED:
    state = CS_SLES_PC_DIVERGED;
    break;
  case CS_SLES_BREAKDOWN:
    state = CS_SLES_PC_BREAKDOWN;
    break;
  case CS_SLES_CONVERGED:
    state = CS_SLES_PC_CONVERGED;
    break;
  default:
    state = CS_SLES_PC_MAX_ITERATION;
  }

  return state;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the MUMPS structure in case of double-precision
 *        computation. Default settings.
 *
 * \param[in]      verbosity   level of verbosity requested
 * \param[in, out] dmumps      pointer to DMUMPS_STRUCT_C
 */
/*----------------------------------------------------------------------------*/

static void
_init_dmumps_settings(int                 verbosity,
                      DMUMPS_STRUC_C     *dmumps)
{
  dmumps->ICNTL(1) = 6;      /* Error output: default value */

  if (verbosity <= 0) {

    dmumps->ICNTL(2) = -1;   /* Rank statistics: default value */
    dmumps->ICNTL(3) = -1;   /* No global information printed */
    dmumps->ICNTL(4) = 1;    /* Only error message printed */
    dmumps->ICNTL(11) = 0;   /* No error analysis */

  }
  else {

    if (verbosity == 1) {
      dmumps->ICNTL(2) = -1;    /* Rank statistics: default value */
      dmumps->ICNTL(3) = 6;     /* Global information: default value */
      dmumps->ICNTL(4) = 1;     /* Only error messages printed */
    }
    else if (verbosity == 2) {

      dmumps->ICNTL(2) = 6;    /* Rank statistics: default value */
      dmumps->ICNTL(3) = 6;    /* Global information: default value */
      dmumps->ICNTL(4) = 2;    /* Verbosity level: default value */

    }
    else { /* verbosity > 2 */

      dmumps->ICNTL(2) = 6;    /* Rank statistics: default value */
      dmumps->ICNTL(3) = 6;    /* Global information: default value */
      dmumps->ICNTL(4) = 4;    /* All messages are printed */
      dmumps->ICNTL(11) = 2;   /* Main error analysis */

    }

  }

  dmumps->ICNTL(5) = 0;    /* 0: assembled / 1: elemental */
  dmumps->ICNTL(20) = 0;   /* 0: dense RHS on rank 0 */
  dmumps->ICNTL(21) = 0;   /* 0: dense solution array on rank 0 */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; MSR matrix as input; no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_msr_dmumps(int                   verbosity,
            const cs_matrix_t    *a,
            DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym == 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  dmumps->n = (MUMPS_INT)n_rows;

  /* Count number of entries (filtering zero or nearly zero values).
   * No modification for the diagonal entries. */

  dmumps->nnz = (MUMPS_INT8)(n_rows);

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++)
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
      if (fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold)
        dmumps->nnz += 1;

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, double);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    dmumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->a[row_id] = (double)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  double  *_a = dmumps->a + n_rows;
  cs_lnum_t  count = 0;

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {
      assert(a_col_ids[i] < n_rows);

      if (fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold) {

        _irn[count] = row_num;
        _jcn[count] = (MUMPS_INT)(a_col_ids[i] + 1);
        _a[count] = (double)x_val[i];
        count++;

      }

    } /* Loop on columns */

  } /* Loop on rows */

  assert(count + n_rows == dmumps->nnz);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system in case of parallel computation.
 *        Case of double-precision structure; MSR matrix as input; no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_parall_msr_dmumps(int                   verbosity,
                   const cs_matrix_t    *a,
                   DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym == 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 3;   /* 3 = distributed matrix is given */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /* Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t  *halo = cs_matrix_get_halo(a);
  const cs_gnum_t  *row_g_id = cs_matrix_get_block_row_g_id(a);

  bool  have_perio = _have_perio(halo);
  CS_UNUSED(have_perio);

  cs_gnum_t  n_g_rows = n_rows;
  cs_parall_counter(&n_g_rows, 1);
  dmumps->n = n_g_rows; /* Global number of rows */

  /* Count number of entries (filtering zero or nearly zero values).
   * No modification for the diagonal entries. */

  dmumps->nnz_loc = (MUMPS_INT8)(n_rows);

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++)
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
      if (fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold)
        dmumps->nnz_loc += 1;

  /* Allocate local arrays */

  BFT_MALLOC(dmumps->irn_loc, dmumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn_loc, dmumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(dmumps->a_loc, dmumps->nnz_loc, double);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    dmumps->irn_loc[row_id] = (MUMPS_INT)row_gnum;
    dmumps->jcn_loc[row_id] = (MUMPS_INT)row_gnum;
    dmumps->a_loc[row_id] = (double)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn_loc + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn_loc + n_rows;
  double  *_a = dmumps->a_loc + n_rows;
  cs_lnum_t  count = 0;

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      if (fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold) {

        _irn[count] = (MUMPS_INT)row_gnum;
        _jcn[count] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
        _a[count] = (double)x_val[i];
        count++;

      }

    } /* Loop on columns */

  } /* Loop on rows */

  assert(count + n_rows == dmumps->nnz_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_native_dmumps(int                   verbosity,
               const cs_matrix_t    *a,
               DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym == 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  /* Retrieve local arrays associated to the current matrix */

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  assert(symmetric == false);
  dmumps->n = (MUMPS_INT)n_rows;
  dmumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, double);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    dmumps->irn[i] = (MUMPS_INT)(i + 1);
    dmumps->jcn[i] = (MUMPS_INT)(i + 1);
    dmumps->a[i] = (double)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  double  *_a = dmumps->a + n_rows;

  cs_lnum_t  count = 0;
  for (cs_lnum_t i = 0; i < n_faces; i++) {

    MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
    MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

    if (c0_id < n_rows) {
      _irn[count] = c0_id + 1;
      _jcn[count] = (MUMPS_INT)(c1_id + 1);
      _a[count] = (double)x_val[2*i];
      count++;
    }

    if (c1_id < n_rows) {
      _irn[count] = c1_id + 1;
      _jcn[count] = (MUMPS_INT)(c0_id + 1);
      _a[count] = (double)x_val[2*i+1];
      count++;
    }

  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system in case of parallel computation.
 *        Case of double-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_parall_native_dmumps(int                   verbosity,
                      const cs_matrix_t    *a,
                      DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym == 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 3;   /* 3 = distributed matrix is given */

  /* Retrieve the matrix arrays */

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  assert(symmetric == false);

  /* Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t  *halo = cs_matrix_get_halo(a);
  const cs_gnum_t  *row_g_id = cs_matrix_get_block_row_g_id(a);

  bool  have_perio = _have_perio(halo);
  CS_UNUSED(have_perio);

  cs_gnum_t  n_g_rows = n_rows;
  cs_parall_counter(&n_g_rows, 1);
  dmumps->n = n_g_rows;

  dmumps->nnz_loc = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(dmumps->irn_loc, dmumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn_loc, dmumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(dmumps->a_loc, dmumps->nnz_loc, double);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    cs_gnum_t  row_gnum = row_g_id[i] + 1;
    dmumps->irn_loc[i] = (MUMPS_INT)(row_gnum);
    dmumps->jcn_loc[i] = (MUMPS_INT)(row_gnum);
    dmumps->a_loc[i] = (double)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn_loc + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn_loc + n_rows;
  double  *_a = dmumps->a_loc + n_rows;

  cs_lnum_t  count = 0;
  for (cs_lnum_t i = 0; i < n_faces; i++) {

    MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
    MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

    if (c0_id < n_rows) {
      _irn[count] = row_g_id[c0_id] + 1;
      _jcn[count] = (MUMPS_INT)(row_g_id[c1_id] + 1);
      _a[count] = (double)x_val[2*i];
      count++;
    }

    if (c1_id < n_rows) {
      _irn[count] = row_g_id[c1_id] + 1;
      _jcn[count] = (MUMPS_INT)(row_g_id[c0_id] + 1);
      _a[count] = (double)x_val[2*i+1];
      count++;
    }

  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; MSR matrix as input; symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_msr_sym_dmumps(int                   verbosity,
                const cs_matrix_t    *a,
                DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym > 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  dmumps->n = (MUMPS_INT)n_rows;
  if (cs_matrix_is_symmetric(a)) /* storage is already symmetric */
    dmumps->nnz = n_rows + a_row_idx[n_rows];

  else {

    /* Count number of entries (filtering zero or nearly zero values).
     * No modification for the diagonal entries. */

    dmumps->nnz = n_rows;

    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++)
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
        if (a_col_ids[i] < row_id &&
            fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold)
          dmumps->nnz += 1;

  }

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, double);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    dmumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->a[row_id] = (double)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  double  *_a = dmumps->a + n_rows;

  if (cs_matrix_is_symmetric(a)) { /* storage is already symmetric */

    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        assert(a_col_ids[i] < n_rows);
        _irn[i] = row_num;
        _jcn[i] = a_col_ids[i] + 1;
        _a[i] = (double)x_val[i];

      } /* Loop on columns */

    } /* Loop on rows */

  }
  else { /* Keep only the lower triangular block */

    cs_lnum_t  count = 0;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        assert(a_col_ids[i] < n_rows);
        MUMPS_INT  col_num = a_col_ids[i] + 1;
        if (col_num < row_num &&
            fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold) {
          _irn[count] = row_num;
          _jcn[count] = col_num;
          _a[count] = (double)x_val[i];
          count++;
        }

      } /* Loop on columns */

    } /* Loop on rows */

    assert(count + n_rows == dmumps->nnz);

  } /* not a symmetric storage */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system in case of parallel computation.
 *        Case of double-precision structure; MSR matrix as input; symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_parall_msr_sym_dmumps(int                   verbosity,
                       const cs_matrix_t    *a,
                       DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym > 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 3;   /* 3 = distributed matrix is given */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t  *halo = cs_matrix_get_halo(a);
  const cs_gnum_t  *row_g_id = cs_matrix_get_block_row_g_id(a);

  bool  have_perio = _have_perio(halo);
  CS_UNUSED(have_perio);        /* TODO */

  cs_gnum_t  n_g_rows = n_rows;
  cs_parall_counter(&n_g_rows, 1);
  dmumps->n = n_g_rows;  /* Global number of rows */

  /* Count number of entries (filtering zero or nearly zero values).
   * No modification for the diagonal entries. */

  if (cs_matrix_is_symmetric(a)) { /* storage is already symmetric */

    dmumps->nnz_loc = n_rows + a_row_idx[n_rows];

  }
  else {

    cs_lnum_t  count = 0;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
        if (row_g_id[a_col_ids[i]] + 1 < row_gnum)
          if (fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold)
            count++;

    } /* Loop on rows */

    dmumps->nnz_loc = n_rows + count;

  }

  /* Allocate local arrays */

  BFT_MALLOC(dmumps->irn_loc, dmumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn_loc, dmumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(dmumps->a_loc, dmumps->nnz_loc, double);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    dmumps->irn_loc[row_id] = (MUMPS_INT)row_gnum;
    dmumps->jcn_loc[row_id] = (MUMPS_INT)row_gnum;
    dmumps->a_loc[row_id] = (double)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn_loc + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn_loc + n_rows;
  double  *_a = dmumps->a_loc + n_rows;

  if (cs_matrix_is_symmetric(a)) { /* storage is already symmetric */

    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        _irn[i] = (MUMPS_INT)row_gnum;
        _jcn[i] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
        _a[i] = (double)x_val[i];

      } /* Loop on columns */

    } /* Loop on rows */

  }
  else { /* Keep only the lower triangular block */

    cs_lnum_t  count = 0;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        const cs_gnum_t  col_gnum = row_g_id[a_col_ids[i]] + 1;
        if (col_gnum < row_gnum &&
            fabs(x_val[i]) > cs_sles_mumps_zero_dthreshold) {

          _irn[count] = (MUMPS_INT)row_gnum;
          _jcn[count] = (MUMPS_INT)col_gnum;
          _a[count] = (double)x_val[i];
          count++;

        }

      } /* Loop on columns */

    } /* Loop on rows */

    assert(count + n_rows == dmumps->nnz_loc);

  } /* Not a symmetric storage */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_native_sym_dmumps(int                   verbosity,
                   const cs_matrix_t    *a,
                   DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym > 0);

  /* Settings */

  _init_dmumps_settings(verbosity, dmumps); /* default settings */

  dmumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  /* Retrieve local arrays associated to the current matrix */

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  dmumps->n = (MUMPS_INT)n_rows;
  dmumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, double);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    dmumps->irn[i] = (MUMPS_INT)(i + 1);
    dmumps->jcn[i] = (MUMPS_INT)(i + 1);
    dmumps->a[i] = (double)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  double  *_a = dmumps->a + n_rows;

  if (symmetric) {

    cs_lnum_t  count = 0;
    for (cs_lnum_t i = 0; i < n_faces; i++) {

      MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
      MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

      if (c0_id < c1_id) {
        _irn[count] = c0_id + 1;
        _jcn[count] = (MUMPS_INT)(c1_id + 1);
        _a[count] = (double)x_val[i];
        count++;
      }
      else {
        assert(c0_id > c1_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (double)x_val[i];
        count++;
      }

    } /* Loop on rows */

  }
  else { /* Native is stored in a non-symmetric way */

    cs_lnum_t  count = 0;
    for (cs_lnum_t i = 0; i < n_faces; i++) {

      MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
      MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

      if (c0_id < c1_id) {
        _irn[count] = c0_id + 1;
        _jcn[count] = (MUMPS_INT)(c1_id + 1);
        _a[count] = (double)x_val[2*i];
        count++;
      }
      else {
        assert(c1_id < c0_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (double)x_val[2*i+1];
        count++;
      }

    } /* Loop on rows */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the MUMPS structure in case of single-precision
 *        computation. Default settings.
 *
 * \param[in]      verbosity   level of verbosity requested
 * \param[in, out] smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_init_smumps_settings(int                 verbosity,
                      SMUMPS_STRUC_C     *smumps)
{
  smumps->ICNTL(1) = 6;      /* Error output: default value */

  if (verbosity <= 0) {

    smumps->ICNTL(2) = -1;   /* Rank statistics: default value */
    smumps->ICNTL(3) = -1;   /* No global information printed */
    smumps->ICNTL(4) = 1;    /* Only error message printed */
    smumps->ICNTL(11) = 0;   /* No error analysis */

  }
  else {

    if (verbosity == 1) {
      smumps->ICNTL(2) = -1;    /* Rank statistics: default value */
      smumps->ICNTL(3) = 6;     /* Global information: default value */
      smumps->ICNTL(4) = 1;     /* Only error messages printed */
    }
    else if (verbosity == 2) {

      smumps->ICNTL(2) = 6;    /* Rank statistics: default value */
      smumps->ICNTL(3) = 6;    /* Global information: default value */
      smumps->ICNTL(4) = 2;    /* Verbosity level: default value */

    }
    else { /* verbosity > 2 */

      smumps->ICNTL(2) = 6;    /* Rank statistics: default value */
      smumps->ICNTL(3) = 6;    /* Global information: default value */
      smumps->ICNTL(4) = 4;    /* All messages are printed */
      smumps->ICNTL(11) = 2;   /* Main error analysis */

    }

  }

  smumps->ICNTL(5) = 0;    /* 0: assembled / 1: elemental */
  smumps->ICNTL(20) = 0;   /* 0: dense RHS on rank 0 */
  smumps->ICNTL(21) = 0;   /* 0: dense solution array on rank 0 */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; MSR matrix as input; no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_msr_smumps(int                   verbosity,
            const cs_matrix_t    *a,
            SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym == 0);

  /* Settings */

  _init_smumps_settings(verbosity, smumps); /* default settings */

  smumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows);

  /* Count number of entries (filtering zero or nearly zero values).
   * No modification for the diagonal entries. */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++)
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
      if (fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold)
        smumps->nnz += 1;

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, float);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    smumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->a[row_id] = (float)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  float  *_a = smumps->a + n_rows;
  cs_lnum_t  count = 0;

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {
      assert(a_col_ids[i] < n_rows);

      if (fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold) {

        _irn[count] = row_num;
        _jcn[count] = (MUMPS_INT)(a_col_ids[i] + 1);
        _a[count] = (float)x_val[i];
        count++;

      }

    } /* Loop on columns */

  } /* Loop on rows */

  assert(count + n_rows == smumps->nnz);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system in case of parallel computation.
 *        Case of single-precision structure; MSR matrix as input; no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_parall_msr_smumps(int                   verbosity,
                   const cs_matrix_t    *a,
                   SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym == 0);

  /* Settings */

  _init_smumps_settings(verbosity, smumps); /* default settings */

  smumps->ICNTL(18) = 3;   /* 3 = distributed matrix is given */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /* Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t  *halo = cs_matrix_get_halo(a);
  const cs_gnum_t  *row_g_id = cs_matrix_get_block_row_g_id(a);

  bool  have_perio = _have_perio(halo);
  CS_UNUSED(have_perio);

  cs_gnum_t  n_g_rows = n_rows;
  cs_parall_counter(&n_g_rows, 1);
  smumps->n = (MUMPS_INT)n_g_rows; /* Global number of rows */

  /* Count number of entries (filtering zero or nearly zero values).
   * No modification for the diagonal entries. */

  smumps->nnz_loc = (MUMPS_INT8)(n_rows);

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++)
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
      if (fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold)
        smumps->nnz_loc += 1;

  /* Allocate local arrays */

  BFT_MALLOC(smumps->irn_loc, smumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(smumps->jcn_loc, smumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(smumps->a_loc, smumps->nnz_loc, float);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    smumps->irn_loc[row_id] = (MUMPS_INT)row_gnum;
    smumps->jcn_loc[row_id] = (MUMPS_INT)row_gnum;
    smumps->a_loc[row_id] = (float)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn_loc + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn_loc + n_rows;
  float  *_a = smumps->a_loc + n_rows;
  cs_lnum_t  count = 0;

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      if (fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold) {

        _irn[count] = (MUMPS_INT)row_gnum;
        _jcn[count] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
        _a[count] = (float)x_val[i];
        count++;

      }

    } /* Loop on columns */

  } /* Loop on rows */

  assert(count + n_rows == smumps->nnz_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_native_smumps(int                   verbosity,
               const cs_matrix_t    *a,
               SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym == 0);

  /* Settings */

  _init_smumps_settings(verbosity, smumps); /* default settings */

  smumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  /* Retrieve local arrays associated to the current matrix */

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  assert(symmetric == false);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, float);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    smumps->irn[i] = (MUMPS_INT)(i + 1);
    smumps->jcn[i] = (MUMPS_INT)(i + 1);
    smumps->a[i] = (float)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  float  *_a = smumps->a + n_rows;

  cs_lnum_t  count = 0;
  for (cs_lnum_t i = 0; i < n_faces; i++) {

    MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
    MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

    if (c0_id < n_rows) {
      _irn[count] = c0_id + 1;
      _jcn[count] = (MUMPS_INT)(c1_id + 1);
      _a[count] = (float)x_val[2*i];
      count++;
    }

    if (c1_id < n_rows) {
      _irn[count] = c1_id + 1;
      _jcn[count] = (MUMPS_INT)(c0_id + 1);
      _a[count] = (float)x_val[2*i+1];
      count++;
    }

  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; MSR matrix as input; symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_msr_sym_smumps(int                   verbosity,
                const cs_matrix_t    *a,
                SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym > 0);

  /* Settings */

  _init_smumps_settings(verbosity, smumps); /* default settings */

  smumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  smumps->n = (MUMPS_INT)n_rows;
  if (cs_matrix_is_symmetric(a)) /* storage is already symmetric */
    smumps->nnz = n_rows + a_row_idx[n_rows];

  else {

    smumps->nnz = (MUMPS_INT8)(n_rows);

    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++)
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
        if (a_col_ids[i] < row_id &&
            fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold)
          smumps->nnz += 1;

  }

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, float);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    smumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->a[row_id] = (float)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  float  *_a = smumps->a + n_rows;

  if (cs_matrix_is_symmetric(a)) { /* storage is already symmetric */

    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        assert(a_col_ids[i] < n_rows);
        _irn[i] = row_num;
        _jcn[i] = a_col_ids[i] + 1;
        _a[i] = (float)x_val[i];

      } /* Loop on columns */

    } /* Loop on rows */

  }
  else { /* Keep only the lower triangular block */

    cs_lnum_t  count = 0;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        assert(a_col_ids[i] < n_rows);
        MUMPS_INT  col_num = a_col_ids[i] + 1;
        if (col_num < row_num &&
            fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold) {
          _irn[count] = row_num;
          _jcn[count] = col_num;
          _a[count] = (float)x_val[i];
          count++;
        }

      } /* Loop on columns */

    } /* Loop on rows */

    assert(count + n_rows == smumps->nnz);

  } /* not a symmetric storage */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system in case of parallel computation.
 *        Case of single-precision structure; MSR matrix as input; symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_parall_msr_sym_smumps(int                   verbosity,
                       const cs_matrix_t    *a,
                       SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym > 0);

  /* Settings */

  _init_smumps_settings(verbosity, smumps); /* default settings */

  smumps->ICNTL(18) = 3;   /* 3 = distributed matrix is given */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t  *halo = cs_matrix_get_halo(a);
  const cs_gnum_t  *row_g_id = cs_matrix_get_block_row_g_id(a);

  bool  have_perio = _have_perio(halo);
  CS_UNUSED(have_perio);        /* TODO */

  cs_gnum_t  n_g_rows = n_rows;
  cs_parall_counter(&n_g_rows, 1);
  smumps->n = n_g_rows;  /* Global number of rows */

  /* Count number of entries (filtering zero or nearly zero values).
   * No modification for the diagonal entries. */

  if (cs_matrix_is_symmetric(a)) { /* storage is already symmetric */

    smumps->nnz_loc = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

  }
  else {

    cs_lnum_t  count = 0;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++)
        if (row_g_id[a_col_ids[i]] + 1 < row_gnum)
          if (fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold)
            count++;

    } /* Loop on rows */

    smumps->nnz_loc = (MUMPS_INT8)(n_rows + count);

  }

  /* Allocate local arrays */

  BFT_MALLOC(smumps->irn_loc, smumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(smumps->jcn_loc, smumps->nnz_loc, MUMPS_INT);
  BFT_MALLOC(smumps->a_loc, smumps->nnz_loc, float);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    smumps->irn_loc[row_id] = (MUMPS_INT)row_gnum;
    smumps->jcn_loc[row_id] = (MUMPS_INT)row_gnum;
    smumps->a_loc[row_id] = (float)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn_loc + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn_loc + n_rows;
  float  *_a = smumps->a_loc + n_rows;

  if (cs_matrix_is_symmetric(a)) { /* storage is already symmetric */

    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        _irn[i] = (MUMPS_INT)row_gnum;
        _jcn[i] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
        _a[i] = (float)x_val[i];

      } /* Loop on columns */

    } /* Loop on rows */

  }
  else { /* Keep only the lower triangular block */

    cs_lnum_t  count = 0;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        const cs_gnum_t  col_gnum = row_g_id[a_col_ids[i]] + 1;
        if (col_gnum < row_gnum &&
            fabs(x_val[i]) > cs_sles_mumps_zero_fthreshold) {

          _irn[count] = (MUMPS_INT)row_gnum;
          _jcn[count] = (MUMPS_INT)col_gnum;
          _a[count] = (float)x_val[i];
          count++;

        }

      } /* Loop on columns */

    } /* Loop on rows */

    assert(count + n_rows == smumps->nnz_loc);

  } /* Not a symmetric storage */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; Native matrix as input; symmetry
 *
 * \param[in]       verbosity   level of verbosity
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_native_sym_smumps(int                   verbosity,
                   const cs_matrix_t    *a,
                   SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym > 0);

  /* Settings */

  _init_smumps_settings(verbosity, smumps); /* default settings */

  smumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 (sequential run) */

  /* Retrieve local arrays associated to the current matrix */

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, float);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    smumps->irn[i] = (MUMPS_INT)(i + 1);
    smumps->jcn[i] = (MUMPS_INT)(i + 1);
    smumps->a[i] = (float)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  float  *_a = smumps->a + n_rows;

  if (symmetric) {

    cs_lnum_t  count = 0;
    for (cs_lnum_t i = 0; i < n_faces; i++) {

      MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
      MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

      if (c0_id < c1_id) {
        _irn[count] = c0_id + 1;
        _jcn[count] = (MUMPS_INT)(c1_id + 1);
        _a[count] = (float)x_val[i];
        count++;
      }
      else {
        assert(c0_id > c1_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (float)x_val[i];
        count++;
      }

    } /* Loop on rows */

    assert(count == n_faces);
  }
  else { /* Native is stored in a non-symmetric way */

    cs_lnum_t  count = 0;
    for (cs_lnum_t i = 0; i < n_faces; i++) {

      MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
      MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

      if (c0_id < c1_id) {
        _irn[count] = c0_id + 1;
        _jcn[count] = (MUMPS_INT)(c1_id + 1);
        _a[count] = (float)x_val[2*i];
        count++;
      }
      else {
        assert(c1_id < c0_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (float)x_val[2*i+1];
        count++;
      }

    } /* Loop on rows */

    assert(count == n_faces);
  }

  /* Update the nnz */

  smumps->nnz = (MUMPS_INT8)(n_rows + n_faces);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if one has to do a new analysis/facto due to an error rising
 *        from the MUMPS feedback. Handle the lack of memory.
 *        Case of DMUMPS structure
 *
 * \param[in, out] c     structure to manage MUMPS execution
 */
/*----------------------------------------------------------------------------*/

static bool
_try_again_dmumps(cs_sles_mumps_t     *c)
{
  DMUMPS_STRUC_C  *mumps = c->mumps_struct;

  int  infog1 = mumps->INFOG(1);

  if (infog1 ==  -8 || infog1 ==  -9 || infog1 == -14 ||
      infog1 == -15 || infog1 == -17 || infog1 == -20) {

    mumps->ICNTL(14) *= 2;  /* Double the portion of memory used in some
                               parts of MUMPS */

    c->n_tries += 1;

    if (c->n_tries > cs_sles_mumps_n_max_tries)
      return false;

    else {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s: DMUMPS has detected a lack of memory.\n"
                    " A new analysis/factorization cycle is going to start.",
                    __func__);

      return true;

    }

  }
  else if (infog1 == -13 || infog1 == -19) {

    mumps->ICNTL(23) = 0;

    c->n_tries += 1;

    if (c->n_tries > cs_sles_mumps_n_max_tries)
      return false;

    else {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s: DMUMPS has detected a lack of memory.\n"
                    " A new analysis/factorization cycle is going to start.",
                    __func__);

      return true;

    }

  }

  return false; /* No need to perform a new computation */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic settings derived from the high-level interface between
 *        code_saturne and MUMPS. These settings are done before the analysis
 *        step and they can still be modified by the user thanks to the function
 *        \ref cs_user_sles_mumps_hook
 *        Case of double-precision MUMPS
 *
 * \param[in]      type    type of factorization to handle
 * \param[in]      slesp   pointer to the related cs_param_sles_t structure
 * \param[in, out] mumps   pointer to a DMUMPS_STRUC_C struct.
 */
/*----------------------------------------------------------------------------*/

static void
_automatic_dmumps_settings_before_analysis(cs_sles_mumps_type_t     type,
                                           const cs_param_sles_t   *slesp,
                                           DMUMPS_STRUC_C          *mumps)
{
  CS_NO_WARN_IF_UNUSED(type);

  cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

  /* Set the algorithm for the analysis step: renumbering and graph
     manipulations */

  switch (mumpsp->analysis_algo) {

  case CS_PARAM_SLES_ANALYSIS_AMD:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 0;
    break;

  case CS_PARAM_SLES_ANALYSIS_QAMD:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 6;
    break;

  case CS_PARAM_SLES_ANALYSIS_PORD:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 4;
    break;

  case CS_PARAM_SLES_ANALYSIS_SCOTCH:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 3;
    mumps->ICNTL(58) = 2;  /* Acceleration of the symbolic factorization */
    break;

  case CS_PARAM_SLES_ANALYSIS_PTSCOTCH:
    mumps->ICNTL(28) = 2;  /* parallel analysis */
    mumps->ICNTL(29) = 1;
    mumps->ICNTL(58) = 0;   /* No symbolic factorization */
    break;

  case CS_PARAM_SLES_ANALYSIS_METIS:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 5;
    break;

  case CS_PARAM_SLES_ANALYSIS_PARMETIS:
    mumps->ICNTL(28) = 2;  /* parallel analysis */
    mumps->ICNTL(29) = 2;
    mumps->ICNTL(58) = 2;  /* Acceleration of the symbolic factorization */
    break;

  default: /* CS_PARAM_SLES_ANALYSIS_AUTO: */
    mumps->ICNTL(7) = 7;
    break;

  } /* Type of algorithm for the analysis step */

  /* Analysis by block if requested */

  if (mumpsp->block_analysis > 1)
    mumps->ICNTL(15) = -mumpsp->block_analysis;

  /* More advanced optimized settings */

  if (mumpsp->advanced_optim) {

    mumps->KEEP(268) = -2;  /* Relaxed pivoting for large enough frontal
                               matrices */

    if (cs_glob_n_ranks > 1 || cs_glob_n_threads > 1) {

      mumps->KEEP(370) = 1; /* Better memory consumption prediction */
      mumps->KEEP(371) = 1; /* Advanced optimization */

    }

    if (cs_glob_n_threads > 1)
      mumps->KEEP(401) = 1; /* Activate openMP tree parallelism (L0-threads) */

  }

  /* Iterative refinement */

  if (mumpsp->ir_steps > 0)
    mumps->ICNTL(10) = -mumpsp->ir_steps; /* Fixed number of iterations */

  /* BLR compression */

  if (fabs(mumpsp->blr_threshold) > FLT_MIN) {

    mumps->ICNTL(35) = 2; /* Activate BLR algo (facto + solve) */
    mumps->CNTL(7) = fabs(mumpsp->blr_threshold); /* Compression rate */

    if (mumpsp->blr_threshold < 0)
      mumps->ICNTL(36) = 0; /* Variante de BLR0 */
    else
      mumps->ICNTL(36) = 1; /* Variante de BLR1 */

    if (mumpsp->mem_usage == CS_PARAM_SLES_MEMORY_CONSTRAINED) {

      mumps->ICNTL(37) = 1; /* Memory compression but time consumming */
      mumps->ICNTL(40) = 1; /* Memory compression mixed precision */
      mumps->ICNTL(49) = 1; /* Memory cleaning of the workspace */

    }
    else {

      mumps->ICNTL(37) = 0; /* No memory compression to save time */
      mumps->ICNTL(40) = 0; /* No memory compression to save time */

    }

  }

  /* Memory workspace */

  if (mumpsp->mem_coef > 0)
    mumps->ICNTL(14) = mumpsp->mem_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic settings derived from the high-level interface between
 *        code_saturne and MUMPS. These settings are done before the
 *        factorization step and they can still be modified by the user thanks
 *        to the function \ref cs_user_sles_mumps_hook
 *        Case of double-precision MUMPS
 *
 * \param[in]      slesp   pointer to the related cs_param_sles_t structure
 * \param[in, out] mumps   pointer to a DMUMPS_STRUC_C struct.
 */
/*----------------------------------------------------------------------------*/

static void
_automatic_dmumps_settings_before_facto(const cs_param_sles_t   *slesp,
                                        DMUMPS_STRUC_C          *mumps)
{
  cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

  if (mumpsp->mem_usage == CS_PARAM_SLES_MEMORY_CPU_DRIVEN) {

    unsigned long  max_estimated_mem = mumps->INFOG(16); /* in MB */
    if (mumps->ICNTL(35) > 1) /* BLR activated */
      max_estimated_mem = mumps->INFOG(36);

    /* In-core usage only up to now.
     * Compute the free memory in RAM */

    struct sysinfo  sys;
    short  status = sysinfo(&sys);

    unsigned long  max_mem_space = max_estimated_mem;
    if (status >= 0)
      max_mem_space = sys.freeram/(1024*1024*sys.mem_unit);

    unsigned long  mem_space = 90*CS_MIN(2*max_estimated_mem,
                                         max_mem_space) / 100;

    mumps->ICNTL(23) = mem_space;

    if (slesp->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT, " MUMPS:"
                    " Estimation of the memory requirement: %lu --> %lu MB\n",
                    max_estimated_mem, mem_space);

  }
  else {

    unsigned long  max_estimated_mem = mumps->INFOG(16); /* in MB */
    if (mumps->ICNTL(35) > 1) /* BLR activated */
      max_estimated_mem = mumps->INFOG(36);

    if (slesp->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT, " MUMPS:"
                    " Estimation of the memory requirement: %lu MB\n",
                    max_estimated_mem);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if one has to do a new analysis/facto due to an error rising
 *        from the MUMPS feedback. Handle the lack of memory.
 *        Case of a SMUMPS structure
 *
 * \param[in, out] c     structure to manage MUMPS execution
 */
/*----------------------------------------------------------------------------*/

static bool
_try_again_smumps(cs_sles_mumps_t     *c)
{
  SMUMPS_STRUC_C  *mumps = c->mumps_struct;

  int  infog1 = mumps->INFOG(1);

  if (infog1 ==  -8 || infog1 ==  -9 || infog1 == -14 ||
      infog1 == -15 || infog1 == -17 || infog1 == -20) {

    mumps->ICNTL(14) *= 2;  /* Double the portion of memory used in some
                               parts of MUMPS */

    c->n_tries += 1;

    if (c->n_tries > cs_sles_mumps_n_max_tries)
      return false;

    else {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s: SMUMPS has detected a lack of memory.\n"
                    " A new analysis/factorization cycle is going to start.",
                    __func__);

      return true;

    }

  }
  else if (infog1 == -13 || infog1 == -19) {

    mumps->ICNTL(23) = 0;

    c->n_tries += 1;

    if (c->n_tries > cs_sles_mumps_n_max_tries)
      return false;

    else {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s: SMUMPS has detected a lack of memory.\n"
                    " A new analysis/factorization cycle is going to start.",
                    __func__);

      return true;

    }

  }

  return false; /* No need to perform a new computation */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic settings derived from the high-level interface between
 *        code_saturne and MUMPS. These settings are done before the analysis
 *        step and they can still be modified by the user thanks to the function
 *        \ref cs_user_sles_mumps_hook
 *        Case of single-precision MUMPS
 *
 * \param[in]      type    type of factorization to handle
 * \param[in]      slesp   pointer to the related cs_param_sles_t structure
 * \param[in, out] mumps   pointer to a SMUMPS_STRUC_C struct.
 */
/*----------------------------------------------------------------------------*/

static void
_automatic_smumps_settings_before_analysis(cs_sles_mumps_type_t     type,
                                           const cs_param_sles_t   *slesp,
                                           SMUMPS_STRUC_C          *mumps)
{
  CS_NO_WARN_IF_UNUSED(type);

  cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

  /* Set the algorithm for the analysis step: renumbering and graph
     manipulations */

  switch (mumpsp->analysis_algo) {

  case CS_PARAM_SLES_ANALYSIS_AMD:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 0;
    break;

  case CS_PARAM_SLES_ANALYSIS_QAMD:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 6;
    break;

  case CS_PARAM_SLES_ANALYSIS_PORD:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 4;
    break;

  case CS_PARAM_SLES_ANALYSIS_SCOTCH:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 3;
    mumps->ICNTL(58) = 2;  /* Acceleration of the symbolic factorization */
    break;

  case CS_PARAM_SLES_ANALYSIS_PTSCOTCH:
    mumps->ICNTL(28) = 2;  /* parallel analysis */
    mumps->ICNTL(29) = 1;
    mumps->ICNTL(58) = 0;  /* No symbolic factorization */
    break;

  case CS_PARAM_SLES_ANALYSIS_METIS:
    mumps->ICNTL(28) = 1;  /* sequential analysis */
    mumps->ICNTL(7) = 5;
    break;

  case CS_PARAM_SLES_ANALYSIS_PARMETIS:
    mumps->ICNTL(28) = 2;  /* parallel analysis */
    mumps->ICNTL(29) = 2;
    mumps->ICNTL(58) = 2;  /* Acceleration of the symbolic factorization */
    break;

  default: /* CS_PARAM_SLES_ANALYSIS_AUTO: */
    mumps->ICNTL(7) = 7;
    break;

  } /* Type of algorithm for the analysis step */

  /* Analysis by block if requested */

  if (mumpsp->block_analysis > 1)
    mumps->ICNTL(15) = -mumpsp->block_analysis;

  /* More advanced optimized settings */

  if (mumpsp->advanced_optim) {

    mumps->KEEP(268) = -2;  /* relaxed pivoting for large enough frontal
                               matrices */

    if (cs_glob_n_ranks > 1 || cs_glob_n_threads > 1) {

      mumps->KEEP(370) = 1; /* Better memory consumption prediction */
      mumps->KEEP(371) = 1; /* Advanced optimization */

    }

    if (cs_glob_n_threads > 1)
      mumps->KEEP(401) = 1; /* Activate openMP tree parallelism (L0-threads) */

  }

  /* Iterative refinement */

  if (mumpsp->ir_steps > 0)
    mumps->ICNTL(10) = -mumpsp->ir_steps; /* Fixed number of iterations */

  /* BLR compression */

  if (fabs(mumpsp->blr_threshold) > FLT_MIN) {

    mumps->ICNTL(35) = 2; /* Activate BLR algo (facto + solve) */
    mumps->CNTL(7) = fabs(mumpsp->blr_threshold); /* Compression rate */

    if (mumpsp->blr_threshold < 0)
      mumps->ICNTL(36) = 0; /* Variante de BLR0 */
    else
      mumps->ICNTL(36) = 1; /* Variante de BLR1 */

    if (mumpsp->mem_usage == CS_PARAM_SLES_MEMORY_CONSTRAINED) {

      mumps->ICNTL(37) = 1; /* Memory compression but time consumming */
      mumps->ICNTL(40) = 1; /* Memory compression mixed precision */
      mumps->ICNTL(49) = 1; /* Memory cleaning of the workspace */

    }
    else {

      mumps->ICNTL(37) = 0; /* No memory compression to save time */
      mumps->ICNTL(40) = 0; /* No memory compression to save time */

    }

  }

  /* Memory workspace */

  if (mumpsp->mem_coef > 0)
    mumps->ICNTL(14) = mumpsp->mem_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic settings derived from the high-level interface between
 *        code_saturne and MUMPS. These settings are done before the
 *        factorization step and they can still be modified by the user thanks
 *        to the function \ref cs_user_sles_mumps_hook
 *        Case of single-precision MUMPS
 *
 * \param[in]      slesp   pointer to the related cs_param_sles_t structure
 * \param[in, out] mumps   pointer to a SMUMPS_STRUC_C struct.
 */
/*----------------------------------------------------------------------------*/

static void
_automatic_smumps_settings_before_facto(const cs_param_sles_t   *slesp,
                                        SMUMPS_STRUC_C          *mumps)
{
  cs_param_sles_mumps_t  *mumpsp = slesp->context_param;

  if (mumpsp->mem_usage == CS_PARAM_SLES_MEMORY_CPU_DRIVEN) {

    unsigned long  max_estimated_mem = mumps->INFOG(16); /* in MB */
    if (mumps->ICNTL(35) > 1) /* BLR activated */
      max_estimated_mem = mumps->INFOG(36);

    /* In-core usage only up to now.
     * Compute the free memory in RAM */

    struct sysinfo  sys;
    short  status = sysinfo(&sys);

    unsigned long  max_mem_space = max_estimated_mem;
    if (status >= 0)
      max_mem_space = sys.freeram/(1024*1024*sys.mem_unit);

    unsigned long  mem_space = 90*CS_MIN(2*max_estimated_mem,
                                         max_mem_space) / 100;

    mumps->ICNTL(23) = mem_space;

    if (slesp->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT, " MUMPS:"
                    " Estimation of the memory requirement: %lu --> %lu MB\n",
                    max_estimated_mem, mem_space);

  }
  else {

    unsigned long  max_estimated_mem = mumps->INFOG(16); /* in MB */
    if (mumps->ICNTL(35) > 1) /* BLR activated */
      max_estimated_mem = mumps->INFOG(36);

    if (slesp->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT, " MUMPS:"
                    " Estimation of the memory requirement: %lu MB\n",
                    max_estimated_mem);

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for advanced user settings of a MUMPS solver.
 *        This function is called two times during the setup stage.
 *        1. Before the analysis step
 *        2. Before the factorization step
 *
 * One can recover the MUMPS step through the "job" member.
 * MUMPS_JOB_ANALYSIS or MUMPS_JOB_FACTORIZATION
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that structure should
 * not be temporary (i.e. local);
 *
 * \param[in]      slesp    pointer to the related cs_param_sles_t structure
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] pmumps   pointer to DMUMPS_STRUC_C or SMUMPS_STRUC_C struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_sles_mumps_hook(const cs_param_sles_t   *slesp,
                        void                    *context,
                        void                    *pmumps)
{
  CS_UNUSED(slesp);
  CS_UNUSED(context);
  CS_UNUSED(pmumps);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate a MUMPS linear system solver for a given field
 *        or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_mumps_create.
 *
 * Note that this function returns a pointer directly to the sparse direct
 * solver management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or NULL
 * \param[in]      slesp         pointer to a cs_param_sles_t structure
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to newly created sparse direct solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_mumps_t *
cs_sles_mumps_define(int                            f_id,
                     const char                    *name,
                     const cs_param_sles_t         *slesp,
                     cs_sles_mumps_setup_hook_t    *setup_hook,
                     void                          *context)
{
  cs_sles_mumps_t * c = cs_sles_mumps_create(slesp, setup_hook, context);

  cs_sles_define(f_id,
                 name,
                 c,
                 "cs_sles_mumps_t",
                 cs_sles_mumps_setup,
                 cs_sles_mumps_solve,
                 cs_sles_mumps_free,
                 cs_sles_mumps_log,
                 cs_sles_mumps_copy,
                 cs_sles_mumps_destroy);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a preconditioner structure relying on MUMPS solver
 *
 * \param[in]      slesp         pointer to a cs_param_sles_t structure
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_mumps_pc_create(const cs_param_sles_t       *slesp)
{
  if (slesp == NULL)
    return NULL;

  assert(slesp->precond == CS_PARAM_PRECOND_MUMPS);

  cs_sles_mumps_t  *c = cs_sles_mumps_create(slesp,
                                             cs_user_sles_mumps_hook,
                                             NULL);
  assert(c->is_pc == true);

  cs_sles_pc_t *pc = cs_sles_pc_define(c,
                                       _mumps_pc_get_type,
                                       _mumps_pc_setup,
                                       NULL, /* tolerance */
                                       _mumps_pc_apply,
                                       cs_sles_mumps_free,
                                       cs_sles_mumps_log,
                                       cs_sles_mumps_copy,
                                       cs_sles_mumps_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create MUMPS linear system solver info and context.
 *
 * \param[in]      slesp         pointer to a cs_param_sles_t structure
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to associated linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_mumps_t *
cs_sles_mumps_create(const cs_param_sles_t       *slesp,
                     cs_sles_mumps_setup_hook_t  *setup_hook,
                     void                        *context)
{
  cs_sles_mumps_t  *c = NULL;

  _n_mumps_systems += 1;

  BFT_MALLOC(c, 1, cs_sles_mumps_t);

  c->type = _set_type(slesp);

  c->n_tries = 0;
  c->n_setups = 0;
  c->n_solves = 0;

  CS_TIMER_COUNTER_INIT(c->t_setup);
  CS_TIMER_COUNTER_INIT(c->t_solve);

  /* Options */

  c->is_pc = _set_pc_usage(slesp);
  c->matrix = NULL;
  c->sles_param = slesp;
  c->hook_context = context;
  c->setup_hook = setup_hook;

  /* Setup data structure */

  c->mumps_struct = NULL;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create MUMPS linear system solver info and context based on existing
 *        info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                      (actual type: cs_sles_mumps_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_mumps_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_mumps_copy(const void   *context)
{
  cs_sles_mumps_t  *d = NULL;

  if (context != NULL) {
    const cs_sles_mumps_t *c = context;
    d = cs_sles_mumps_create(c->sles_param,
                             c->setup_hook,
                             c->hook_context);
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free MUMPS linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to sparse direct solver info and context
 *                           (actual type: cs_sles_mumps_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_mumps_free(void  *context)
{
  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_mumps_t  *c  = context;

  if (c == NULL)
    return;

  if (c->mumps_struct != NULL) {

    if (_is_dmumps(c)) {

      DMUMPS_STRUC_C  *dmumps = c->mumps_struct;

      dmumps->job = MUMPS_JOB_END;
      dmumps_c(dmumps);

      if (cs_glob_n_ranks == 1) {

        BFT_FREE(dmumps->irn);
        BFT_FREE(dmumps->jcn);
        BFT_FREE(dmumps->a);

      }
      else {

        BFT_FREE(dmumps->irn_loc);
        BFT_FREE(dmumps->jcn_loc);
        BFT_FREE(dmumps->a_loc);

      }

      BFT_FREE(dmumps);

    }
    else {

      SMUMPS_STRUC_C  *smumps = c->mumps_struct;

      smumps->job = MUMPS_JOB_END;
      smumps_c(smumps);

      if (cs_glob_n_ranks == 1) {

        BFT_FREE(smumps->irn);
        BFT_FREE(smumps->jcn);
        BFT_FREE(smumps->a);

      }
      else {

        BFT_FREE(smumps->irn_loc);
        BFT_FREE(smumps->jcn_loc);
        BFT_FREE(smumps->a_loc);

      }

      BFT_FREE(smumps);

    }

    c->mumps_struct = NULL;

  } /* MUMPS structure is not freed */

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy MUMPS linear system solver info and context.
 *
 * \param[in, out]  context  pointer to sparse direct solver info and context
 *                           (actual type: cs_sles_mumps_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_mumps_destroy(void   **context)
{
  cs_sles_mumps_t *c = (cs_sles_mumps_t *)(*context);
  if (c != NULL) {

    /* Free structure */

    cs_sles_mumps_free(c);
    BFT_FREE(c);
    *context = c;

    _n_mumps_systems -= 1;

  } /* c != NULL */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup MUMPS linear equation solver.
 *
 * \param[in, out]  context    pointer to sparse direct solver info and context
 *                             (actual type: cs_sles_mumps_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_mumps_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity)
{
  CS_UNUSED(name);

  cs_timer_t t0;
  t0 = cs_timer_time();

  /* Sanity checks */

  assert(a != NULL);
  if (cs_matrix_get_diag_block_size(a) > 1 ||
      cs_matrix_get_extra_diag_block_size(a) > 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid matrix structure for MUMPS. No block requested.\n",
              __func__);

  /* Manage the MPI communicator */

  if (cs_glob_n_ranks == 1) { /* sequential run */

#if defined(HAVE_MPI)
    int  flag = 0;
    MPI_Initialized(&flag);

    if (!flag) {
#if   (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
      int mpi_threads;
      MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &mpi_threads);
#else
      MPI_Init(NULL, NULL);
#endif
    }

    if (cs_glob_mpi_comm == MPI_COMM_NULL)
      cs_glob_mpi_comm = MPI_COMM_WORLD;

#endif /* HAVE_MPI */
  }

  /* Begin the setup */

  cs_sles_mumps_t  *c = context;

  c->mumps_struct = NULL;

  /* 1. Initialize the MUMPS structure */
  /* --------------------------------- */

  if (_is_dmumps(c)) {

    /* Sanity checks: DMUMPS_COMPLEX = DMUMPS_REAL = double
     * (see mumps_c_types.h) */

    assert(sizeof(double) == sizeof(DMUMPS_COMPLEX));
    assert(sizeof(double) == sizeof(DMUMPS_REAL));

    DMUMPS_STRUC_C  *dmumps = NULL;
    BFT_MALLOC(dmumps, 1, DMUMPS_STRUC_C);

    dmumps->job = MUMPS_JOB_INIT;
    dmumps->par = 1;      /* all ranks are working */

    if (c->type == CS_SLES_MUMPS_DOUBLE_LU)
      dmumps->sym = 0;
    else if (c->type == CS_SLES_MUMPS_DOUBLE_LDLT_SPD)
      dmumps->sym = 1;
    else if (c->type == CS_SLES_MUMPS_DOUBLE_LDLT_SYM)
      dmumps->sym = 2;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid type of MUMPS settings for double-precision.\n",
                __func__);

#if defined(HAVE_MPI)
    dmumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else
    /* Not used in this case and set to the default value given by the MUMPS
       documentation */

    dmumps->comm_fortran = USE_COMM_WORLD;
#endif

    dmumps_c(dmumps); /* first call to MUMPS: Initialization */

    /* Set the MUMPS pointer */

    c->mumps_struct = dmumps;

  }
  else {

    /* Sanity checks: SMUMPS_COMPLEX = SMUMPS_REAL = float
     * (see mumps_c_types.h) */

    assert(sizeof(float) == sizeof(SMUMPS_COMPLEX));
    assert(sizeof(float) == sizeof(SMUMPS_REAL));

    SMUMPS_STRUC_C  *smumps = NULL;
    BFT_MALLOC(smumps, 1, SMUMPS_STRUC_C);

    smumps->job = MUMPS_JOB_INIT;
    smumps->par = 1;       /* all ranks are working */
    smumps->sym = 0;

    if (c->type == CS_SLES_MUMPS_SINGLE_LU)
      smumps->sym = 0;
    else if (c->type == CS_SLES_MUMPS_SINGLE_LDLT_SPD)
      smumps->sym = 1;
    else if (c->type == CS_SLES_MUMPS_SINGLE_LDLT_SYM)
      smumps->sym = 2;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid type of MUMPS settings for single-precision.\n",
                __func__);


#if defined(HAVE_MPI)
    smumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else

    /* Not used in this case and set to the default value given by the MUMPS
       documentation */

    smumps->comm_fortran = USE_COMM_WORLD;
#endif

    smumps_c(smumps); /* first call to MUMPS: Initialization */

    /* Set the MUMPS pointer */

    c->mumps_struct = smumps;

  }

  /* 2. Fill the MUMPS structure before the analysis step */
  /* ---------------------------------------------------- */

  const cs_matrix_type_t  cs_mat_type = cs_matrix_get_type(a);

  switch (c->type) {

  case CS_SLES_MUMPS_DOUBLE_LU:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_dmumps(verbosity, a, c->mumps_struct);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _parall_native_dmumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_dmumps(verbosity, a, c->mumps_struct);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_dmumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_SLES_MUMPS_DOUBLE_LDLT_SPD:
  case CS_SLES_MUMPS_DOUBLE_LDLT_SYM:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_sym_dmumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_sym_dmumps(verbosity, a, c->mumps_struct);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_sym_dmumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_SLES_MUMPS_SINGLE_LU:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_smumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_smumps(verbosity, a, c->mumps_struct);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_smumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_SLES_MUMPS_SINGLE_LDLT_SPD:
  case CS_SLES_MUMPS_SINGLE_LDLT_SYM:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_sym_smumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_sym_smumps(verbosity, a, c->mumps_struct);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_sym_smumps(verbosity, a, c->mumps_struct);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: MUMPS is not set as a solver.\n"
              " Please check your settings.", __func__);

  } /* End of switch */

  /* 3. Analysis and factorization */
  /* ----------------------------- */

  MUMPS_INT  infog1, infog2;

  c->n_tries = 0;    /* reset the number of tries */

  if (_is_dmumps(c)) {

    DMUMPS_STRUC_C  *dmumps = c->mumps_struct;

    do {

      /* Analysis step */
      /* ------------- */

      dmumps->job = MUMPS_JOB_ANALYSIS;

      _automatic_dmumps_settings_before_analysis(c->type, c->sles_param,
                                                 dmumps);

      /* Window to enable advanced user settings (before analysis) */

      if (c->setup_hook != NULL)
        c->setup_hook(c->sles_param, c->hook_context, dmumps);

      dmumps_c(dmumps);

      /* Factorization step */
      /* ------------------ */

      dmumps->job = MUMPS_JOB_FACTORIZATION;

      _automatic_dmumps_settings_before_facto(c->sles_param, dmumps);

      /* Window to enable advanced user settings (before factorization) */

      if (c->setup_hook != NULL)
        c->setup_hook(c->sles_param, c->hook_context, dmumps);

      dmumps_c(dmumps);

      /* Feedback */

      infog1 = dmumps->INFOG(1);
      infog2 = dmumps->INFOG(2);

    } while (_try_again_dmumps(c));

  }
  else {

    SMUMPS_STRUC_C  *smumps = c->mumps_struct;

    do {

      /* Analysis step */
      /* ------------- */

      smumps->job = MUMPS_JOB_ANALYSIS;

      _automatic_smumps_settings_before_analysis(c->type, c->sles_param,
                                                 smumps);

      /* Window to enable advanced user settings (before analysis) */

      if (c->setup_hook != NULL)
        c->setup_hook(c->sles_param, c->hook_context, smumps);

      smumps_c(smumps);

      /* Factorization step */
      /* ------------------ */

      smumps->job = MUMPS_JOB_FACTORIZATION;

      _automatic_smumps_settings_before_facto(c->sles_param, smumps);

      /* Window to enable advanced user settings (before factorization) */

      if (c->setup_hook != NULL)
        c->setup_hook(c->sles_param, c->hook_context, smumps);

      smumps_c(smumps);

      /* Feedback */

      infog1 = smumps->INFOG(1);
      infog2 = smumps->INFOG(2);

    } while(_try_again_smumps(c));

  } /* single-precision case */

  /* Check the feedback given by MUMPS */

  if (cs_glob_rank_id < 1) {

    if (infog1 < 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n MUMPS feedback error code: INFOG(1)=%d, INFOG(2)=%d\n",
                    infog1, infog2);
      bft_error(__FILE__, __LINE__, 0,
                "%s: Error detected during the analysis/factorization step.\n"
                "%s: INFOG(1)=%d; INFOG(2)=%d\n"
                " Please refer to the MUMPS documentation to get a more"
                " detailed feedback.\n",
                __func__, __func__, infog1, infog2);
    }
    else {
      if (verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n MUMPS feedback code: INFOG(1)=%d, INFOG(2)=%d\n",
                      infog1, infog2);
    }

  } /* rank_id = 0 */

  /* Update returned values */

  c->n_setups += 1;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call MUMPS linear equation solver.
 *
 * \param[in, out]  context        pointer to sparse direct solver info and
 *                                 context (actual type: cs_sles_mumps_t  *)
 * \param[in]       name           pointer to system name
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if NULL)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_mumps_solve(void                *context,
                    const char          *name,
                    const cs_matrix_t   *a,
                    int                  verbosity,
                    double               precision,
                    double               r_norm,
                    int                 *n_iter,
                    double              *residual,
                    const cs_real_t     *rhs,
                    cs_real_t           *vx,
                    size_t               aux_size,
                    void                *aux_vectors)
{
  CS_UNUSED(precision);
  CS_UNUSED(r_norm);
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_sles_mumps_t  *c = context;

  if (c->mumps_struct == NULL)
    cs_sles_mumps_setup(c, name, a, verbosity);

  MUMPS_INT  infog1 = 0;
  cs_timer_t t0;
  t0 = cs_timer_time();

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  cs_fp_exception_disable_trap();

  if (_is_dmumps(c)) {

    /* MUMPS with double-precision arrays */
    /* ---------------------------------- */

    DMUMPS_STRUC_C  *dmumps = c->mumps_struct;
    assert(dmumps != NULL);

    /* 1. Build the RHS */

    if (cs_glob_n_ranks == 1) { /* Sequential run */

      assert(n_rows == dmumps->n);
      dmumps->nrhs = 1;
      cs_array_real_copy(n_rows, rhs, vx);
      dmumps->rhs = vx;

    }
    else { /* Parallel computation */

      assert(cs_glob_n_ranks > 1);

      /* Gather on the rank 0 (= host rank for MUMPS) the global RHS array */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = dmumps->n;

      double  *glob_rhs = NULL;
      if (cs_glob_rank_id == root_rank)
        BFT_MALLOC(glob_rhs, n_g_rows, double);

      cs_parall_gather_r(root_rank, n_rows, n_g_rows, rhs, glob_rhs);

      dmumps->rhs = glob_rhs;

    }

    /* 2. Resolution */

    dmumps->job = MUMPS_JOB_SOLVE;
    dmumps_c(dmumps);
    infog1 = dmumps->INFOG(1);     /* feedback */
    *residual = dmumps->RINFOG(11); /* scaled residual */

    /* 3. Post-resolution operations */

    if (cs_glob_n_ranks == 1)
      dmumps->rhs = NULL;

    else {

      /* Scatter operation (solution is stored in the RHS array.
       * Elements in glob_rhs belonging to a distant rank are sent back to
       * this rank
       */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = dmumps->n;
      double  *glob_rhs = dmumps->rhs;

      cs_parall_scatter_r(root_rank, n_rows, n_g_rows, glob_rhs, vx);

      if (cs_glob_rank_id == root_rank)
        BFT_FREE(glob_rhs);
      dmumps->rhs = NULL;

    }

  }
  else {

    /* MUMPS with single-precision arrays */
    /* ---------------------------------- */

    SMUMPS_STRUC_C  *smumps = c->mumps_struct;
    assert(smumps != NULL);

    /* 1. Build the RHS */

    if (cs_glob_n_ranks == 1) { /* Sequential run */

      assert(n_rows == smumps->n);
      smumps->nrhs = 1;
      BFT_MALLOC(smumps->rhs, n_rows, float);

      /* The MUMPS structure stores the RHS with the type SMUMPS_COMPLEX */

      for (cs_lnum_t i = 0; i < n_rows; i++)
        smumps->rhs[i] = (float)rhs[i];

    }
    else {

      assert(cs_glob_n_ranks > 1);

      /* Gather on the rank 0 (= host rank for MUMPS) the global RHS array */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = smumps->n;
      float  *glob_rhs = NULL;
      if (cs_glob_rank_id == root_rank)
        BFT_MALLOC(glob_rhs, n_g_rows, float);

      /* Use vx (an initial guess is not useful for a direct solver) to define
       * a single-precision rhs
       */

      float  *_svx = (float *)vx;

#     pragma omp parallel for if (n_rows > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_rows; i++)
        _svx[i] = (float)rhs[i];

      cs_parall_gather_f(root_rank, n_rows, n_g_rows, _svx, glob_rhs);

      smumps->rhs = glob_rhs;

    }

    /* 2. Resolution */

    smumps->job = MUMPS_JOB_SOLVE;
    smumps_c(smumps);
    infog1 = smumps->INFOG(1);     /* feedback */
    *residual = smumps->RINFOG(11); /* scaled residual */

    /* 3. Post-resolution operations */

    if (cs_glob_n_ranks == 1) {

      /* Solution is stored inside rhs. Copy/cast into vx */

      for (cs_lnum_t i = 0; i < n_rows; i++)
        vx[i] = (cs_real_t)smumps->rhs[i];

      BFT_FREE(smumps->rhs);

    }
    else {

      /* Scatter operation (solution is stored in the RHS array).
       * Elements in glob_rhs belonging to a distant rank are sent back to
       * this rank.
       */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = smumps->n;
      float  *glob_rhs = smumps->rhs;
      float  *_svx = (float *)vx;

      cs_parall_scatter_f(root_rank, n_rows, n_g_rows, glob_rhs, _svx);

      /* avoid overwritting since sizeof(float) should lower than
         sizeof(cs_real_t) */
      for (cs_lnum_t i = n_rows-1; i > -1; i--)
        vx[i] = (cs_real_t)_svx[i];

      BFT_FREE(glob_rhs);

    }

    smumps->rhs = NULL;

  } /* Single or double-precision algorithm */

  cs_fp_exception_restore_trap();

  /* Output */

  cs_sles_convergence_state_t cvg = CS_SLES_CONVERGED;
  if (infog1 < 0) {

    cvg = CS_SLES_BREAKDOWN;

    if (verbosity > 0) {
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT, "%s: MUMPS Feedback: INFOG(1)=%d\n",
                    __func__, infog1);
    }

  }

  *n_iter = 1;
  c->n_solves += 1;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to sparse direct solver info and context
 *                       (actual type: cs_sles_mumps_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_mumps_log(const void  *context,
                  cs_log_t     log_type)
{
  const cs_sles_mumps_t  *c = context;

  char sym_type_name[32];
  char storage_type_name[32];

  switch(c->type) {
  case CS_SLES_MUMPS_DOUBLE_LU:
    strncpy(sym_type_name, "non-symmetric", 31);
    strncpy(storage_type_name, "double-precision", 31);
    break;
  case CS_SLES_MUMPS_DOUBLE_LDLT_SPD:
    strncpy(sym_type_name, "symmetric; SPD", 31);
    strncpy(storage_type_name, "double-precision", 31);
    break;
  case CS_SLES_MUMPS_DOUBLE_LDLT_SYM:
    strncpy(sym_type_name, "general symmetric", 31);
    strncpy(storage_type_name, "double-precision", 31);
    break;
  case CS_SLES_MUMPS_SINGLE_LU:
    strncpy(sym_type_name, "non-symmetric", 31);
    strncpy(storage_type_name, "single-precision", 31);
    break;
  case CS_SLES_MUMPS_SINGLE_LDLT_SPD:
    strncpy(sym_type_name, "symmetric; SPD", 31);
    strncpy(storage_type_name, "single-precision", 31);
    break;
  case CS_SLES_MUMPS_SINGLE_LDLT_SYM:
    strncpy(sym_type_name, "general symmetric", 31);
    strncpy(storage_type_name, "single-precision", 31);
    break;

  default:
    strncpy(sym_type_name, "unknown", 31);
    strncpy(storage_type_name, "unknown", 31);
  }

  sym_type_name[31] = '\0';
  storage_type_name[31] = '\0';

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  "  Solver type:                       MUMPS %s\n"
                  "    Storage:                           %s\n"
                  "    Symm type:                         %s\n",
                  MUMPS_VERSION, storage_type_name, sym_type_name);

  }
  else if (log_type == CS_LOG_PERFORMANCE) {

    if (c->is_pc)
      cs_log_printf(log_type,
                    _("\n"
                      "  Preconditioner type:           MUMPS\n"
                      "  Number of setups:              %12d\n"
                      "  Number of solves:              %12d\n"
                      "  Total setup time:              %12.3f\n"
                      "  Total solution time:           %12.3f\n"),
                    c->n_setups, c->n_solves,
                    c->t_setup.nsec*1e-9, c->t_solve.nsec*1e-9);
    else
      cs_log_printf(log_type,
                    _("\n"
                      "  Solver type:                   MUMPS\n"
                      "  Number of setups:              %12d\n"
                      "  Number of solves:              %12d\n"
                      "  Total setup time:              %12.3f\n"
                      "  Total solution time:           %12.3f\n"),
                    c->n_setups, c->n_solves,
                    c->t_setup.nsec*1e-9, c->t_solve.nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on MUMPS library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_mumps_library_info(cs_log_t  log_type)
{
  cs_log_printf(log_type, "    MUMPS %s\n", MUMPS_VERSION);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
