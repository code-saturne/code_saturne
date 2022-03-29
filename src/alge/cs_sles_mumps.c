/*============================================================================
 * Sparse Linear Equation Solvers using MUMPS (a sparse direct solver library)
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

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

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_sles_mumps_setup_t {

  /* single precision structure */

  SMUMPS_STRUC_C      *smumps;

  /* double precision structure */

  DMUMPS_STRUC_C      *dmumps;

} cs_sles_mumps_setup_t;

struct _cs_sles_mumps_t {

  /* Performance data */

  int                  n_setups;      /* Number of times system setup */
  int                  n_solves;      /* Number of times system solved since
                                       * it's a direct solver:
                                       * n_solves = n_iterations_tot */

  cs_timer_counter_t   t_setup;       /* Total setup (factorization) */
  cs_timer_counter_t   t_solve;       /* Total time used */

  /* Additional setup options */

  const cs_param_sles_t       *sles_param;    /* set of parameter for SLES */

  void                        *hook_context;  /* Optional user context */
  cs_sles_mumps_setup_hook_t  *setup_hook;    /* Post setup function */

  /* Setup data context */

  cs_sles_mumps_setup_t       *setup_data;

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static int  _n_mumps_systems = 0;

/*============================================================================
 * Static inline private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the MUMPS structure in case of double-precision
 *        computation.
 *
 * \param[in]      verbosity   level of verbosity requested
 * \param[in, out] dmumps      pointer to DMUMPS_STRUCT_C
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

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
  dmumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

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

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {
      assert(a_col_ids[i] < n_rows);

      _irn[i] = row_num;
      _jcn[i] = (MUMPS_INT)(a_col_ids[i] + 1);
      _a[i] = (double)x_val[i];

    } /* Loop on columns */

  } /* Loop on rows */

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

  /* Allocate local arrays */

  dmumps->nnz_loc = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

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

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      _irn[i] = (MUMPS_INT)row_gnum;
      _jcn[i] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
      _a[i] = (double)x_val[i];

    } /* Loop on columns */

  } /* Loop on rows */

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
  dmumps->CNTL(1) = 0.0;   /* No pivoting (quicker) */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  dmumps->n = (MUMPS_INT)n_rows;
  if (cs_matrix_is_symmetric(a)) /* storage is already symmetric */
    dmumps->nnz = n_rows + a_row_idx[n_rows];
  else
    dmumps->nnz = n_rows + a_row_idx[n_rows]/2;

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
        if (col_num < row_num) {
          _irn[count] = row_num;
          _jcn[count] = col_num;
          _a[count] = (double)x_val[i];
          count++;
        }

      } /* Loop on columns */

    } /* Loop on rows */

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
  dmumps->CNTL(1) = 0.0;   /* No pivoting (quicker) */

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

  /* Allocate local arrays */

  if (cs_matrix_is_symmetric(a)) /* storage is already symmetric */
    dmumps->nnz_loc = n_rows + a_row_idx[n_rows];
  else
    dmumps->nnz_loc = n_rows + a_row_idx[n_rows]/2;

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

        if (a_col_ids[i] < row_id) {
          _irn[count] = (MUMPS_INT)row_gnum;
          _jcn[count] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
          _a[count] = (double)x_val[i];
          count++;
        }

      } /* Loop on columns */

    } /* Loop on rows */

  } /* not a symmetric storage */
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
  dmumps->CNTL(1) = 0.0;   /* No pivoting (quicker) */

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
  smumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

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

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {
      assert(a_col_ids[i] < n_rows);

      _irn[i] = row_num;
      _jcn[i] = (MUMPS_INT)(a_col_ids[i] + 1);
      _a[i] = (float)x_val[i];

    } /* Loop on columns */

  } /* Loop on rows */

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

  /* Allocate local arrays */

  smumps->nnz_loc = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

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

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      _irn[i] = (MUMPS_INT)row_gnum;
      _jcn[i] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
      _a[i] = (float)x_val[i];

    } /* Loop on columns */

  } /* Loop on rows */

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
  smumps->CNTL(1) = 0.0;   /* No pivoting (quicker) */

  /* Retrieve local arrays associated to the current matrix */

  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  /*  Fill the MUMPS matrix */

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

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

  cs_lnum_t  count = n_rows;
  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      assert(a_col_ids[i] < n_rows);
      MUMPS_INT  col_num = a_col_ids[i] + 1;
      if (col_num < row_num) {
        _irn[count] = row_num;
        _jcn[count] = col_num;
        _a[count] = (float)x_val[i];
        count++;
      }

    } /* Loop on columns */

  } /* Loop on rows */

  /* TODO: optimization when the storage is already symmetric */
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
  smumps->CNTL(1) = 0.0;   /* No pivoting (quicker) */
  smumps->ICNTL(6) = 0;    /* No column permutation */

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

  /* Allocate local arrays */

  if (cs_matrix_is_symmetric(a)) /* storage is already symmetric */
    smumps->nnz_loc = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);
  else
    smumps->nnz_loc = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]/2);

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

    cs_lnum_t  count = n_rows;
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

      const cs_gnum_t  row_gnum = row_g_id[row_id] + 1;
      for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

        if (a_col_ids[i] < row_id) {
          _irn[count] = (MUMPS_INT)row_gnum;
          _jcn[count] = (MUMPS_INT)(row_g_id[a_col_ids[i]] + 1);
          _a[count] = (float)x_val[i];
          count++;
        }

      } /* Loop on columns */

    } /* Loop on rows */

  } /* not a symmetric storage */
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
  smumps->CNTL(1) = 0.0;   /* No pivoting (quicker) */

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for user settings of a MUMPS solver.
 *        This function is called at the end of the setup stage.
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that structure should
 * not be temporary (i.e. local);
 *
 * \param[in]      slesp      pointer to the related cs_param_sles_t structure
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] dmumps     pointer to DMUMPS_STRUC_C (double-precision)
 * \param[in, out] smumps     pointer to SMUMPS_STRUC_C (single-precision)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_sles_mumps_hook(const cs_param_sles_t   *slesp,
                        void                    *context,
                        DMUMPS_STRUC_C          *dmumps,
                        SMUMPS_STRUC_C          *smumps)
{
  CS_UNUSED(slesp);
  CS_UNUSED(context);
  CS_UNUSED(dmumps);
  CS_UNUSED(smumps);
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
  c->n_setups = 0;
  c->n_solves = 0;

  CS_TIMER_COUNTER_INIT(c->t_setup);
  CS_TIMER_COUNTER_INIT(c->t_solve);

  /* Options */

  c->sles_param = slesp;
  c->hook_context = context;
  c->setup_hook = setup_hook;

  /* Setup data */

  c->setup_data = NULL;

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

  cs_sles_mumps_setup_t *sd = c->setup_data;

  if (sd != NULL) {

    if (sd->dmumps != NULL) {

      sd->dmumps->job = MUMPS_JOB_END;
      dmumps_c(sd->dmumps);

      if (cs_glob_n_ranks == 1) {

        BFT_FREE(sd->dmumps->irn);
        BFT_FREE(sd->dmumps->jcn);
        BFT_FREE(sd->dmumps->a);

      }
      else {

        BFT_FREE(sd->dmumps->irn_loc);
        BFT_FREE(sd->dmumps->jcn_loc);
        BFT_FREE(sd->dmumps->a_loc);

      }

      BFT_FREE(sd->dmumps);

    }

    if (sd->smumps != NULL) {

      sd->smumps->job = MUMPS_JOB_END;
      smumps_c(sd->smumps);

      if (cs_glob_n_ranks == 1) {

        BFT_FREE(sd->smumps->irn);
        BFT_FREE(sd->smumps->jcn);
        BFT_FREE(sd->smumps->a);

      }
      else {

        BFT_FREE(sd->smumps->irn_loc);
        BFT_FREE(sd->smumps->jcn_loc);
        BFT_FREE(sd->smumps->a_loc);

      }

      BFT_FREE(sd->smumps);

    }

    BFT_FREE(c->setup_data);

  }

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
  cs_sles_mumps_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    BFT_MALLOC(c->setup_data, 1, cs_sles_mumps_setup_t);
    sd = c->setup_data;
  }

  /* 1. Initialize the MUMPS structure */
  /* --------------------------------- */

  if (c->sles_param->solver == CS_PARAM_ITSOL_MUMPS ||
      c->sles_param->solver == CS_PARAM_ITSOL_MUMPS_LDLT) {

    /* Sanity checks: DMUMPS_COMPLEX = DMUMPS_REAL = double
     * (see mumps_c_types.h) */

    assert(sizeof(double) == sizeof(DMUMPS_COMPLEX));
    assert(sizeof(double) == sizeof(DMUMPS_REAL));

    sd->smumps = NULL;        /* Not used anymore */

    BFT_MALLOC(sd->dmumps, 1, DMUMPS_STRUC_C);

    sd->dmumps->job = MUMPS_JOB_INIT;
    sd->dmumps->par = 1;      /* all ranks are working */
    sd->dmumps->sym = 0;

    if (c->sles_param->solver == CS_PARAM_ITSOL_MUMPS_LDLT)
      sd->dmumps->sym = 2;

#if defined(HAVE_MPI)
    sd->dmumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else
    /* Not used in this case and set to the default value given by the MUMPS
       documentation */
    sd->dmumps->comm_fortran = USE_COMM_WORLD;
#endif

    dmumps_c(sd->dmumps); /* first call to MUMPS */

  }
  else if (c->sles_param->solver == CS_PARAM_ITSOL_MUMPS_FLOAT ||
           c->sles_param->solver == CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT) {

    /* Sanity checks: SMUMPS_COMPLEX = SMUMPS_REAL = float
     * (see mumps_c_types.h) */

    assert(sizeof(float) == sizeof(SMUMPS_COMPLEX));
    assert(sizeof(float) == sizeof(SMUMPS_REAL));

    sd->dmumps = NULL;        /* Not used anymore */

    BFT_MALLOC(sd->smumps, 1, SMUMPS_STRUC_C);

    sd->smumps->job = MUMPS_JOB_INIT;
    sd->smumps->par = 1;       /* all ranks are working */
    sd->smumps->sym = 0;

    if (c->sles_param->solver == CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT)
      sd->smumps->sym = 2;

#if defined(HAVE_MPI)
    sd->smumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else
    /* Not used in this case and set to the default value given by the MUMPS
       documentation */
    sd->smumps->comm_fortran = USE_COMM_WORLD;
#endif

    smumps_c(sd->smumps); /* first call to MUMPS */

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of solver for the MUMPS library\n", __func__);

  /* 2. Fill the MUMPS structure before the analysis step */
  /* ---------------------------------------------------- */

  const cs_matrix_type_t  cs_mat_type = cs_matrix_get_type(a);

  switch (c->sles_param->solver) {

  case CS_PARAM_ITSOL_MUMPS:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_dmumps(verbosity, a, sd->dmumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _parall_native_dmumps(verbosity, a, sd->dmumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_dmumps(verbosity, a, sd->dmumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_dmumps(verbosity, a, sd->dmumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_PARAM_ITSOL_MUMPS_LDLT:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_sym_dmumps(verbosity, a, sd->dmumps);
      /* else if (cs_mat_type == CS_MATRIX_NATIVE) */
      /*   _parall_native_sym_dmumps(verbosity, a, sd->dmumps); */
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_sym_dmumps(verbosity, a, sd->dmumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_sym_dmumps(verbosity, a, sd->dmumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_PARAM_ITSOL_MUMPS_FLOAT:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_smumps(verbosity, a, sd->smumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_smumps(verbosity, a, sd->smumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_smumps(verbosity, a, sd->smumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
    if (cs_glob_n_ranks > 1) { /* Parallel computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _parall_msr_sym_smumps(verbosity, a, sd->smumps);
      /* else if (cs_mat_type == CS_MATRIX_NATIVE) */
      /*   _parall_native_sym_smumps(a, sd->smumps); */
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format in parallel", __func__);

    }
    else { /* Sequential computation */

      if (cs_mat_type == CS_MATRIX_MSR)
        _msr_sym_smumps(verbosity, a, sd->smumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _native_sym_smumps(verbosity, a, sd->smumps);
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

  if (sd->smumps == NULL) {

    sd->dmumps->job = MUMPS_JOB_ANALYSIS;

    /* Window to enable advanced user settings (before analysis) */
    if (c->setup_hook != NULL)
      c->setup_hook(c->sles_param, c->hook_context, sd->dmumps, sd->smumps);

    dmumps_c(sd->dmumps);

  }
  else {

    assert(sd->smumps != NULL);
    sd->smumps->job = MUMPS_JOB_ANALYSIS;

    /* Window to enable advanced user settings (before analysis) */
    if (c->setup_hook != NULL)
      c->setup_hook(c->sles_param, c->hook_context, sd->dmumps, sd->smumps);

    smumps_c(sd->smumps);

  }

  /* 3. Factorization */
  /* ---------------- */

  MUMPS_INT  infog1, infog2;

  if (sd->smumps == NULL) {

    sd->dmumps->job = MUMPS_JOB_FACTORIZATION;

    /* Window to enable advanced user settings (before factorization) */
    if (c->setup_hook != NULL)
      c->setup_hook(c->sles_param, c->hook_context, sd->dmumps, sd->smumps);

    dmumps_c(sd->dmumps);

    /* Feedback */
    infog1 = sd->dmumps->INFOG(1);
    infog2 = sd->dmumps->INFOG(2);

  }
  else {

    assert(sd->smumps != NULL);
    sd->smumps->job = MUMPS_JOB_FACTORIZATION;

    /* Window to enable advanced user settings (before factorization) */
    if (c->setup_hook != NULL)
      c->setup_hook(c->sles_param, c->hook_context, sd->dmumps, sd->smumps);

    smumps_c(sd->smumps);

    /* Feedback */
    infog1 = sd->smumps->INFOG(1);
    infog2 = sd->smumps->INFOG(2);

  }

  /* Check the feedback given by MUMPS */

  if (cs_glob_rank_id < 1) {

    if (infog1 < 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n MUMPS feedback error code: INFOG(1)=%d, INFOG(2)=%d\n",
                    infog1, infog2);
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected during the analysis/factorization step",
                __func__);
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
 * \param[in]       r_norm         residue normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residue        residue
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
                    double              *residue,
                    const cs_real_t     *rhs,
                    cs_real_t           *vx,
                    size_t               aux_size,
                    void                *aux_vectors)
{
  CS_UNUSED(precision);
  CS_UNUSED(r_norm);
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_timer_t t0;
  t0 = cs_timer_time();

  MUMPS_INT  infog1 = 0;
  cs_sles_mumps_t  *c = context;
  cs_sles_mumps_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    cs_sles_mumps_setup(c, name, a, verbosity);
    sd = c->setup_data;
  }

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  cs_fp_exception_disable_trap();

  switch (c->sles_param->solver) {

    /* MUMPS with double-precision arrays */
    /* ---------------------------------- */

  case CS_PARAM_ITSOL_MUMPS:
  case CS_PARAM_ITSOL_MUMPS_LDLT:

    /* Sanity checks */

    assert(sd->dmumps != NULL);
    assert(sizeof(double) == sizeof(cs_real_t));

    /* Build the RHS */
    if (cs_glob_n_ranks == 1) { /* sequential run */

      assert(n_rows == sd->dmumps->n);
      sd->dmumps->nrhs = 1;
      memcpy(vx, rhs, n_rows*sizeof(double));
      sd->dmumps->rhs = vx;

    }
    else { /* parallel computation */

      assert(cs_glob_n_ranks > 1);

      /* Gather on the rank 0 (= host rank for MUMPS) the global RHS array */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = sd->dmumps->n;

      double  *glob_rhs = NULL;
      if (cs_glob_rank_id == root_rank)
        BFT_MALLOC(glob_rhs, n_g_rows, double);

      cs_parall_gather_r(root_rank, n_rows, n_g_rows, rhs, glob_rhs);

      sd->dmumps->rhs = glob_rhs;

    }

    /* Resolution */

    sd->dmumps->job = MUMPS_JOB_SOLVE;
    dmumps_c(sd->dmumps);
    infog1 = sd->dmumps->INFOG(1);     /* feedback */
    *residue = sd->dmumps->RINFOG(11); /* scaled residual */

    /* Post-resolution operations */

    if (cs_glob_n_ranks == 1)
      sd->dmumps->rhs = NULL;

    else {

      /* Scatter operation (solution is stored in the RHS array.
       * Elements in glob_rhs belonging to a distant rank are sent back to
       * this rank
       */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = sd->dmumps->n;
      double  *glob_rhs = sd->dmumps->rhs;

      cs_parall_scatter_r(root_rank, n_rows, n_g_rows, glob_rhs, vx);

      if (cs_glob_rank_id == root_rank)
        BFT_FREE(glob_rhs);
      sd->dmumps->rhs = NULL;

    }
    break; /* double-precision */

    /* MUMPS with single-precision arrays */
    /* ---------------------------------- */

  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
  case CS_PARAM_ITSOL_MUMPS_FLOAT:

    /* Sanity checks */

    assert(sd->smumps != NULL);

    /* Build the RHS */

    if (cs_glob_n_ranks == 1) { /* sequential run */

      assert(n_rows == sd->smumps->n);
      sd->smumps->nrhs = 1;
      BFT_MALLOC(sd->smumps->rhs, n_rows, float);

      /* The MUMPS structure stores the RHS with the type SMUMPS_COMPLEX */
      for (cs_lnum_t i = 0; i < n_rows; i++)
        sd->smumps->rhs[i] = (float)rhs[i];

    }
    else {

      assert(cs_glob_n_ranks > 1);

      /* Gather on the rank 0 (= host rank for MUMPS) the global RHS array */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = sd->smumps->n;
      float  *glob_rhs = NULL;
      if (cs_glob_rank_id == root_rank)
        BFT_MALLOC(glob_rhs, n_g_rows, float);

      /* Use vx (an initial guess is not useful for a direct solver) to define
       * a single-precision rhs */
      float  *_svx = (float *)vx;

#     pragma omp parallel for if (n_rows > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_rows; i++)
        _svx[i] = (float)rhs[i];

      cs_parall_gather_f(root_rank, n_rows, n_g_rows, _svx,
                         glob_rhs);

      sd->smumps->rhs = glob_rhs;

    }

    /* Resolution */

    sd->smumps->job = MUMPS_JOB_SOLVE;
    smumps_c(sd->smumps);
    infog1 = sd->smumps->INFOG(1);     /* feedback */
    *residue = sd->smumps->RINFOG(11); /* scaled residual */

    /* Post-resolution operations */

    if (cs_glob_n_ranks == 1) {

      /* Solution is stored inside rhs. Copy/cast into vx */

      for (cs_lnum_t i = 0; i < n_rows; i++)
        vx[i] = (cs_real_t)sd->smumps->rhs[i];

      BFT_FREE(sd->smumps->rhs);

    }
    else {

      /* Scatter operation (solution is stored in the RHS array).
       * Elements in glob_rhs belonging to a distant rank are sent back to
       * this rank. */

      int  root_rank = 0;
      MUMPS_INT  n_g_rows = sd->smumps->n;
      float  *glob_rhs = sd->smumps->rhs;
      float  *_svx = (float *)vx;

      cs_parall_scatter_f(root_rank, n_rows, n_g_rows, glob_rhs, _svx);

      /* avoid overwritting since sizeof(float) should lower than
         sizeof(cs_real_t) */
      for (cs_lnum_t i = n_rows-1; i > -1; i--)
        vx[i] = (cs_real_t)_svx[i];

      BFT_FREE(glob_rhs);

    }

    sd->smumps->rhs = NULL;

    break; /* single-precision */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid solver.\n", __func__);

  }

  cs_fp_exception_restore_trap();

  /* Output */

  cs_sles_convergence_state_t cvg = CS_SLES_CONVERGED;
  if (infog1 < 0)
    cvg = CS_SLES_BREAKDOWN;

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

  switch(c->sles_param->solver) {
  case CS_PARAM_ITSOL_MUMPS:
    strncpy(sym_type_name, "non-symmetric", 31);
    strncpy(storage_type_name, "double-precision", 31);
    break;
  case CS_PARAM_ITSOL_MUMPS_LDLT:
    strncpy(sym_type_name, "symmetric", 31);
    strncpy(storage_type_name, "double-precision", 31);
    break;
  case CS_PARAM_ITSOL_MUMPS_FLOAT:
    strncpy(sym_type_name, "non-symmetric", 31);
    strncpy(storage_type_name, "single-precision", 31);
    break;
  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
    strncpy(sym_type_name, "symmetric", 31);
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
                  "  Solver type:                       MUMPS\n"
                  "    Storage:                           %s\n"
                  "    Symm type:                         %s\n"
                  "    Library version:                   %s\n",
                  storage_type_name, sym_type_name, MUMPS_VERSION);

  }
  else if (log_type == CS_LOG_PERFORMANCE) {

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   MUMPS %s\n"
                    "  Number of setups:              %12d\n"
                    "  Number of solves:              %12d\n"
                    "  Total setup time:              %12.3f\n"
                    "  Total solution time:           %12.3f\n"),
                  MUMPS_VERSION, c->n_setups, c->n_solves,
                  c->t_setup.nsec*1e-9, c->t_solve.nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
