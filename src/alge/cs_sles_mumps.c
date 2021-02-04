/*============================================================================
 * Sparse Linear Equation Solvers using MUMPS (a sparse direct solver library)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

  \brief handling of MUMPS-based linear solvers

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
#define MUMPS_JOB_INIT            -1
#define MUMPS_JOB_FACTSYMBOLIC     1
#define MUMPS_JOB_FACTNUMERIC      2
#define MUMPS_JOB_SOLVE            3
#define MUMPS_JOB_ANALYSE_FACTO    4
#define MUMPS_JOB_END             -2

/* Default number for MPI_COMM_WORLD */
#define USE_COMM_WORLD         -987654

/* Set of macros defined in order to match MUMPS documentation (because of the
 * difference between C/FORTRAN programming language
 */
#define ICNTL(I)   icntl[(I)-1]
#define CNTL(I)    cntl[(I)-1]
#define INFOG(I)   infog[(I)-1]
#define INFO(I)    info[(I)-1]
#define RINFOG(I)  rinfog[(I)-1]
#define RINFO(I)   rinfo[(I)-1]

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
 * Private function definitions
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

static void
_init_dmumps(int                 verbosity,
             DMUMPS_STRUC_C     *dmumps)
{
  /* Manage verbosity and output */
  /* --------------------------- */

  dmumps->ICNTL(1) = 6;      /* Error output: default value */

  if (verbosity <= 0) {

    dmumps->ICNTL(2) = -1;   /* Rank statistics: default value */
    dmumps->ICNTL(3) = -1;   /* No global information printed */
    dmumps->ICNTL(4) = 0;    /* No message printed */
    dmumps->ICNTL(11) = 0;   /* No error analysis */

  }
  else {

    dmumps->ICNTL(2) = 0;    /* Rank statistics: default value */
    dmumps->ICNTL(3) = 6;    /* Global information: default value */
    dmumps->ICNTL(11) = 2;   /* Main error analysis */

    if (verbosity == 1)
      dmumps->ICNTL(4) = 1;  /* Only error messages printed */
    else if (verbosity == 2)
      dmumps->ICNTL(4) = 2;  /* Verbosity level: default value */
    else /* verbosity > 2 */
      dmumps->ICNTL(4) = 4;  /* All messages are printed */

  }

  if (cs_glob_n_ranks == 1) { /* sequential run */

    dmumps->ICNTL(5) = 0;    /* 0: assembled / 1: elemental */
    dmumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 */
    dmumps->ICNTL(20) = 0;   /* 0: dense RHS on rank 0 */
    dmumps->ICNTL(21) = 0;   /* 0: dense solution array on rank 0 */

  }

#if defined(HAVE_MPI)
  dmumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else
  dmumps->comm_fortran = USE_COMM_WORLD; /* Not used in this case and set to
                                            the default value given by the
                                            MUMPS documentation */
#endif

  /* First call: Initialization */
  /* -------------------------- */

  dmumps_c(dmumps);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; MSR matrix as input; no symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_msr_dmumps(const cs_matrix_t    *a,
                  DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym == 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  dmumps->n = (MUMPS_INT)n_rows;
  dmumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, DMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    dmumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->a[row_id] = (DMUMPS_REAL)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  DMUMPS_REAL  *_a = dmumps->a + n_rows;

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {
      assert(a_col_ids[i] < n_rows);

      _irn[i] = row_num;
      _jcn[i] = (MUMPS_INT)(a_col_ids[i] + 1);
      _a[i] = (DMUMPS_REAL)x_val[i];

    }
  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_native_dmumps(const cs_matrix_t    *a,
                     DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym == 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

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
  BFT_MALLOC(dmumps->a, dmumps->nnz, DMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    dmumps->irn[i] = (MUMPS_INT)(i + 1);
    dmumps->jcn[i] = (MUMPS_INT)(i + 1);
    dmumps->a[i] = (DMUMPS_REAL)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  DMUMPS_REAL  *_a = dmumps->a + n_rows;

  cs_lnum_t  count = 0;
  for (cs_lnum_t i = 0; i < n_faces; i++) {

    MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
    MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

    if (c0_id < n_rows) {
      _irn[count] = c0_id + 1;
      _jcn[count] = (MUMPS_INT)(c1_id + 1);
      _a[count] = (DMUMPS_REAL)x_val[2*i];
      count++;
    }

    if (c1_id < n_rows) {
      _irn[count] = c1_id + 1;
      _jcn[count] = (MUMPS_INT)(c0_id + 1);
      _a[count] = (DMUMPS_REAL)x_val[2*i+1];
      count++;
    }

  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; MSR matrix as input; symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_msr_sym_dmumps(const cs_matrix_t    *a,
                      DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym > 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  dmumps->n = (MUMPS_INT)n_rows;
  dmumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, DMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    dmumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    dmumps->a[row_id] = (DMUMPS_REAL)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  DMUMPS_REAL  *_a = dmumps->a + n_rows;

  cs_lnum_t  count = n_rows;
  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      assert(a_col_ids[i] < n_rows);
      MUMPS_INT  col_num = a_col_ids[i] + 1;
      if (col_num < row_num) {
        _irn[count] = row_num;
        _jcn[count] = col_num;
        _a[count] = (DMUMPS_REAL)x_val[i];
        count++;
      }

    }
  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of double-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  dmumps      pointer to DMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_native_sym_dmumps(const cs_matrix_t    *a,
                         DMUMPS_STRUC_C       *dmumps)
{
  assert(dmumps->sym > 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  dmumps->n = (MUMPS_INT)n_rows;
  dmumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(dmumps->irn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->jcn, dmumps->nnz, MUMPS_INT);
  BFT_MALLOC(dmumps->a, dmumps->nnz, DMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    dmumps->irn[i] = (MUMPS_INT)(i + 1);
    dmumps->jcn[i] = (MUMPS_INT)(i + 1);
    dmumps->a[i] = (DMUMPS_REAL)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = dmumps->irn + n_rows;
  MUMPS_INT  *_jcn = dmumps->jcn + n_rows;
  DMUMPS_REAL  *_a = dmumps->a + n_rows;

  if (symmetric) {

    cs_lnum_t  count = 0;
    for (cs_lnum_t i = 0; i < n_faces; i++) {

      MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
      MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

      if (c0_id < c1_id) {
        _irn[count] = c0_id + 1;
        _jcn[count] = (MUMPS_INT)(c1_id + 1);
        _a[count] = (DMUMPS_REAL)x_val[i];
        count++;
      }
      else {
        assert(c0_id > c1_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (DMUMPS_REAL)x_val[i];
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
        _a[count] = (DMUMPS_REAL)x_val[2*i];
        count++;
      }
      else {
        assert(c1_id < c0_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (DMUMPS_REAL)x_val[2*i+1];
        count++;
      }

    } /* Loop on rows */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the MUMPS structure in case of single-precision
 *        computation.
 *
 * \param[in]      verbosity   level of verbosity requested
 * \param[in, out] smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_init_smumps(int                 verbosity,
             SMUMPS_STRUC_C     *smumps)
{
  /* Manage verbosity and output */
  /* --------------------------- */

  smumps->ICNTL(1) = 6;      /* Error output: default value */

  if (verbosity <= 0) {

    smumps->ICNTL(2) = -1;   /* Rank statistics: default value */
    smumps->ICNTL(3) = -1;   /* No global information printed */
    smumps->ICNTL(4) = 0;    /* No message printed */
    smumps->ICNTL(11) = 0;   /* No error analysis */

  }
  else {

    smumps->ICNTL(2) = 0;    /* Rank statistics: default value */
    smumps->ICNTL(3) = 6;    /* Global information: default value */
    smumps->ICNTL(11) = 2;   /* Main error analysis */

    if (verbosity == 1)
      smumps->ICNTL(4) = 1;  /* Only error messages printed */
    else if (verbosity == 2)
      smumps->ICNTL(4) = 2;  /* Verbosity level: default value */
    else /* verbosity > 2 */
      smumps->ICNTL(4) = 4;  /* All messages are printed */

  }

  if (cs_glob_n_ranks == 1) { /* sequential run */

    smumps->ICNTL(5) = 0;    /* 0: assembled / 1: elemental */
    smumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 */
    smumps->ICNTL(20) = 0;   /* 0: dense RHS on rank 0 */
    smumps->ICNTL(21) = 0;   /* 0: dense solution array on rank 0 */

  }

#if defined(HAVE_MPI)
  smumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else
  smumps->comm_fortran = USE_COMM_WORLD; /* Not used in this case and set to
                                            the default value given by the
                                            MUMPS documentation */
#endif

  /* First call: Initialization */
  /* -------------------------- */

  smumps_c(smumps);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; MSR matrix as input; no symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_msr_smumps(const cs_matrix_t    *a,
                  SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym == 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, SMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    smumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->a[row_id] = (SMUMPS_REAL)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  SMUMPS_REAL  *_a = smumps->a + n_rows;

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {
      assert(a_col_ids[i] < n_rows);

      _irn[i] = row_num;
      _jcn[i] = (MUMPS_INT)(a_col_ids[i] + 1);
      _a[i] = (SMUMPS_REAL)x_val[i];

    }
  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; Native matrix as input;
 *        no symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_native_smumps(const cs_matrix_t    *a,
                     SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym == 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  assert(symmetric == false);
  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, SMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    smumps->irn[i] = (MUMPS_INT)(i + 1);
    smumps->jcn[i] = (MUMPS_INT)(i + 1);
    smumps->a[i] = (SMUMPS_REAL)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  SMUMPS_REAL  *_a = smumps->a + n_rows;

  cs_lnum_t  count = 0;
  for (cs_lnum_t i = 0; i < n_faces; i++) {

    MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
    MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

    if (c0_id < n_rows) {
      _irn[count] = c0_id + 1;
      _jcn[count] = (MUMPS_INT)(c1_id + 1);
      _a[count] = (SMUMPS_REAL)x_val[2*i];
      count++;
    }

    if (c1_id < n_rows) {
      _irn[count] = c1_id + 1;
      _jcn[count] = (MUMPS_INT)(c0_id + 1);
      _a[count] = (SMUMPS_REAL)x_val[2*i+1];
      count++;
    }

  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; MSR matrix as input; symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_msr_sym_smumps(const cs_matrix_t    *a,
                      SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym > 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t  *a_row_idx, *a_col_ids;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, SMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    smumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
    smumps->a[row_id] = (SMUMPS_REAL)d_val[row_id];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  SMUMPS_REAL  *_a = smumps->a + n_rows;

  cs_lnum_t  count = n_rows;
  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

    MUMPS_INT  row_num = (MUMPS_INT)(row_id + 1);
    for (cs_lnum_t i = a_row_idx[row_id]; i < a_row_idx[row_id+1]; i++) {

      assert(a_col_ids[i] < n_rows);
      MUMPS_INT  col_num = a_col_ids[i] + 1;
      if (col_num < row_num) {
        _irn[count] = row_num;
        _jcn[count] = col_num;
        _a[count] = (SMUMPS_REAL)x_val[i];
        count++;
      }

    }
  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the linear system.
 *        Case of single-precision structure; Native matrix as input; symmetry
 *
 * \param[in]       a           associated matrix
 * \param[in, out]  smumps      pointer to SMUMPS_STRUC_C
 */
/*----------------------------------------------------------------------------*/

static void
_build_native_sym_smumps(const cs_matrix_t    *a,
                         SMUMPS_STRUC_C       *smumps)
{
  assert(smumps->sym > 0);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  bool  symmetric = false;
  cs_lnum_t  n_faces = 0;
  const cs_lnum_2_t  *face_cells;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_native_arrays(a,
                              &symmetric,
                              &n_faces, &face_cells, &d_val, &x_val);

  smumps->n = (MUMPS_INT)n_rows;
  smumps->nnz = (MUMPS_INT8)(n_rows + 2*n_faces);

  BFT_MALLOC(smumps->irn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->jcn, smumps->nnz, MUMPS_INT);
  BFT_MALLOC(smumps->a, smumps->nnz, SMUMPS_REAL);

  /* Add diagonal entries */

  for (cs_lnum_t i = 0; i < n_rows; i++) {

    smumps->irn[i] = (MUMPS_INT)(i + 1);
    smumps->jcn[i] = (MUMPS_INT)(i + 1);
    smumps->a[i] = (SMUMPS_REAL)d_val[i];

  }

  /* Extra-diagonal entries */

  MUMPS_INT  *_irn = smumps->irn + n_rows;
  MUMPS_INT  *_jcn = smumps->jcn + n_rows;
  SMUMPS_REAL  *_a = smumps->a + n_rows;

  if (symmetric) {

    cs_lnum_t  count = 0;
    for (cs_lnum_t i = 0; i < n_faces; i++) {

      MUMPS_INT  c0_id = (MUMPS_INT)(face_cells[i][0]);
      MUMPS_INT  c1_id = (MUMPS_INT)(face_cells[i][1]);

      if (c0_id < c1_id) {
        _irn[count] = c0_id + 1;
        _jcn[count] = (MUMPS_INT)(c1_id + 1);
        _a[count] = (SMUMPS_REAL)x_val[i];
        count++;
      }
      else {
        assert(c0_id > c1_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (SMUMPS_REAL)x_val[i];
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
        _a[count] = (SMUMPS_REAL)x_val[2*i];
        count++;
      }
      else {
        assert(c1_id < c0_id);
        _irn[count] = c1_id + 1;
        _jcn[count] = (MUMPS_INT)(c0_id + 1);
        _a[count] = (SMUMPS_REAL)x_val[2*i+1];
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
      else
        bft_error(__FILE__, __LINE__, 0, "%s: Not yet implemented", __func__);

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
      else
        bft_error(__FILE__, __LINE__, 0, "%s: Not yet implemented", __func__);

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
  if (cs_matrix_get_diag_block_size(a)[0] > 1 ||
      cs_matrix_get_extra_diag_block_size(a)[0] > 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid matrix structure for MUMPS. No block requested.\n",
              __func__);

  const cs_halo_t  *halo = cs_matrix_get_halo(a);

  if (halo != NULL)
    if (halo->n_transforms > 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: No periodicity requested with MUMPS.\n", __func__);

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
  else { /* Parallel computation */

    bft_error(__FILE__, __LINE__, 0,
              " %s: Parallel implementation not available yet.\n",
              __func__);

  }

  /* Begin the setup */

  cs_sles_mumps_t  *c = context;
  cs_sles_mumps_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    BFT_MALLOC(c->setup_data, 1, cs_sles_mumps_setup_t);
    sd = c->setup_data;
  }

  int _verbosity = c->sles_param->verbosity;
  if (_verbosity < 0)
    _verbosity = verbosity;

  const cs_matrix_type_t  cs_mat_type = cs_matrix_get_type(a);

  switch (c->sles_param->solver) {

  case CS_PARAM_ITSOL_MUMPS:
    {
      sd->smumps = NULL;

      BFT_MALLOC(sd->dmumps, 1, DMUMPS_STRUC_C);

      sd->dmumps->job = MUMPS_JOB_INIT;
      sd->dmumps->par = 1;       /* all ranks are working */
      sd->dmumps->sym = 0;

      /* Initialization */

      _init_dmumps(_verbosity, sd->dmumps);

      /* Set the linear system */
      if (cs_mat_type == CS_MATRIX_MSR)
        _build_msr_dmumps(a, sd->dmumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _build_native_dmumps(a, sd->dmumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_PARAM_ITSOL_MUMPS_FLOAT:
    {
      sd->dmumps = NULL;
      BFT_MALLOC(sd->smumps, 1, SMUMPS_STRUC_C);

      sd->smumps->job = MUMPS_JOB_INIT;
      sd->smumps->par = 1;       /* all ranks are working */
      sd->smumps->sym = 0;

      _init_smumps(_verbosity, sd->smumps);

      /* Set the linear system */
      if (cs_mat_type == CS_MATRIX_MSR)
        _build_msr_smumps(a, sd->smumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _build_native_smumps(a, sd->smumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);
    }
    break;

  case CS_PARAM_ITSOL_MUMPS_LDLT:
    {
      sd->smumps = NULL;
      BFT_MALLOC(sd->dmumps, 1, DMUMPS_STRUC_C);

      sd->dmumps->job = MUMPS_JOB_INIT;
      sd->dmumps->par = 1;       /* all ranks are working */
      sd->dmumps->sym = 2;

      _init_dmumps(_verbosity, sd->dmumps);

      /* Set the linear system */
      if (cs_mat_type == CS_MATRIX_MSR)
        _build_msr_sym_dmumps(a, sd->dmumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _build_native_sym_dmumps(a, sd->dmumps);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid matrix format", __func__);

    }
    break;

  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
    {
      sd->dmumps = NULL;
      BFT_MALLOC(sd->smumps, 1, SMUMPS_STRUC_C);

      sd->smumps->job = MUMPS_JOB_INIT;
      sd->smumps->par = 1;       /* all ranks are working */
      sd->smumps->sym = 2;

      _init_smumps(_verbosity, sd->smumps);

      /* Set the linear system */
      if (cs_mat_type == CS_MATRIX_MSR)
        _build_msr_sym_smumps(a, sd->smumps);
      else if (cs_mat_type == CS_MATRIX_NATIVE)
        _build_native_sym_smumps(a, sd->smumps);
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

  /* Window to enable advanced user settings */
  if (c->setup_hook != NULL)
    c->setup_hook(c->sles_param, c->hook_context,
                  sd->dmumps, sd->smumps);

  /* Update return values */

  c->n_setups += 1;

  /* Second call: analysis and factorization */
  /* --------------------------------------- */

  MUMPS_INT  infog1, infog2;

  if (sd->smumps == NULL) {

    sd->dmumps->job = MUMPS_JOB_ANALYSE_FACTO;
    dmumps_c(sd->dmumps);

    /* Feedback */
    infog1 = sd->dmumps->INFOG(1);
    infog2 = sd->dmumps->INFOG(2);

  }
  else {

    assert(sd->smumps != NULL);
    sd->smumps->job = MUMPS_JOB_ANALYSE_FACTO;
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
                " %s: Error detected during the anaylis/factorization step",
                __func__);
    }
    else {
      if (verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n MUMPS feedback code: INFOG(1)=%d, INFOG(2)=%d\n",
                      infog1, infog2);
    }

  } /* rank_id = 0 */

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
 * \param[in]       rotation_mode  halo update option for rotational periodicity
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
                    cs_halo_rotation_t   rotation_mode,
                    double               precision,
                    double               r_norm,
                    int                 *n_iter,
                    double              *residue,
                    const cs_real_t     *rhs,
                    cs_real_t           *vx,
                    size_t               aux_size,
                    void                *aux_vectors)
{
  CS_UNUSED(rotation_mode);
  CS_UNUSED(precision);
  CS_UNUSED(r_norm);
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_timer_t t0;
  t0 = cs_timer_time();

  MUMPS_INT  infog1;
  cs_sles_mumps_t  *c = context;
  cs_sles_mumps_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    cs_sles_mumps_setup(c, name, a, verbosity);
    sd = c->setup_data;
  }

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  cs_fp_exception_disable_trap();

  switch (c->sles_param->solver) {

  case CS_PARAM_ITSOL_MUMPS:
  case CS_PARAM_ITSOL_MUMPS_LDLT:
    assert(sd->dmumps != NULL);
    assert(n_rows == sd->dmumps->n);
    assert(sizeof(cs_real_t) == sizeof(DMUMPS_REAL));

    /* Build the RHS */
    if (cs_glob_n_ranks == 1) { /* sequential run */

      sd->dmumps->nrhs = 1;
      memcpy(vx, rhs, n_rows*sizeof(cs_real_t));
      sd->dmumps->rhs = vx;

    }
    else {

      assert(cs_glob_n_ranks > 1);

      sd->dmumps->nloc_rhs = sd->dmumps->n;
      sd->dmumps->lrhs_loc = sd->dmumps->n;
      BFT_MALLOC(sd->dmumps->rhs_loc, sd->dmumps->n, DMUMPS_REAL);
      BFT_MALLOC(sd->dmumps->irhs_loc, sd->dmumps->n, MUMPS_INT);

      /* TODO */

    }

    /* Resolution */

    sd->dmumps->job = MUMPS_JOB_SOLVE;
    dmumps_c(sd->dmumps);
    infog1 = sd->dmumps->INFOG(1); /* feedback */
    *residue = sd->dmumps->RINFOG(11); /* scaled residual */

    /* Free buffers */
    if (cs_glob_n_ranks == 1)
      sd->dmumps->rhs = NULL;

    else {

      BFT_FREE(sd->dmumps->rhs_loc);
      BFT_FREE(sd->dmumps->irhs_loc);

      sd->dmumps->rhs_loc = NULL;
      sd->dmumps->irhs_loc = NULL;

    }
    break;

  case CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT:
  case CS_PARAM_ITSOL_MUMPS_FLOAT:
    /* Sanity checks */
    assert(sd->smumps != NULL);
    assert(n_rows == sd->smumps->n);

    /* Build the RHS */
    if (cs_glob_n_ranks == 1) { /* sequential run */

      sd->smumps->nrhs = 1;
      BFT_MALLOC(sd->smumps->rhs, n_rows, SMUMPS_REAL);
      for (cs_lnum_t i = 0; i < n_rows; i++)
        sd->smumps->rhs[i] = (float)rhs[i];

    }
    else {

      assert(cs_glob_n_ranks > 1);

      sd->smumps->nloc_rhs = sd->smumps->n;
      sd->smumps->lrhs_loc = sd->smumps->n;
      BFT_MALLOC(sd->smumps->rhs_loc, sd->smumps->n, SMUMPS_REAL);
      BFT_MALLOC(sd->smumps->irhs_loc, sd->smumps->n, MUMPS_INT);

      /* TODO */

    }

    /* Resolution */

    sd->smumps->job = MUMPS_JOB_SOLVE;
    smumps_c(sd->smumps);
    infog1 = sd->smumps->INFOG(1); /* feedback */
    *residue = sd->smumps->RINFOG(11); /* scaled residual */

    /* Free buffers */
    if (cs_glob_n_ranks == 1) {

      /* Solution is stored inside rhs. Copy/cast into vx */
      for (cs_lnum_t i = 0; i < n_rows; i++)
        vx[i] = (cs_real_t)sd->smumps->rhs[i];

      BFT_FREE(sd->smumps->rhs);
      sd->smumps->rhs = NULL;

    }
    else {

      BFT_FREE(sd->smumps->rhs_loc);
      BFT_FREE(sd->smumps->irhs_loc);

      sd->smumps->rhs_loc = NULL;
      sd->smumps->irhs_loc = NULL;

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid solver.\n", __func__);

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
                  "    Symm type:                         %s\n",
                  storage_type_name, sym_type_name);

  }
  else if (log_type == CS_LOG_PERFORMANCE) {

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   MUMPS\n"
                    "  Number of setups:              %12d\n"
                    "  Number of solves:              %12d\n"
                    "  Total setup time:              %12.3f\n"
                    "  Total solution time:           %12.3f\n"),
                  c->n_setups, c->n_solves,
                  c->t_setup.wall_nsec*1e-9, c->t_solve.wall_nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
