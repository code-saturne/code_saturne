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

  DMUMPS_STRUC_C      *mumps;           /* Linear solver context: including
                                         * matrix, rhs and settings */

  void                *cctx;            /* convergence context */

} cs_sles_mumps_setup_t;

struct _cs_sles_mumps_t {

  /* Performance data */

  int                  n_setups;           /* Number of times system setup */
  int                  n_solves;           /* Number of times system solved
                                            * Since it's a direct solver:
                                            * n_solves = n_iterations_tot */

  cs_timer_counter_t   t_setup;            /* Total setup (factorization) */
  cs_timer_counter_t   t_solve;            /* Total time used */

  /*Additional setup options */

  int                          sym;           /* 0: unsymmetric
                                               * 1: S.P.D.
                                               * 2: symmetric */
  int                          verbosity;

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * User function prototypes
 *============================================================================*/
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
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] mumps    pointer to DMUMPS_STRUCT_C structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_sles_mumps_hook(void               *context,
                        DMUMPS_STRUC_C     *mumps)
{
  CS_UNUSED(context);
  CS_UNUSED(mumps);
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
 * \param[in]      sym           type of matrix (unsymmetric, SPD, symmetric)
 * \param[in]      verbosity     level of verbosity
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to newly created sparse direct solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_mumps_t *
cs_sles_mumps_define(int                          f_id,
                     const char                  *name,
                     int                          sym,
                     int                          verbosity,
                     cs_sles_mumps_setup_hook_t  *setup_hook,
                     void                        *context)
{
  cs_sles_mumps_t * c = cs_sles_mumps_create(sym,
                                             verbosity,
                                             setup_hook,
                                             context);

  cs_sles_t *sc = cs_sles_define(f_id,
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
 * \param[in]      sym           type of matrix (unsymmetric, SPD, symmetric)
 * \param[in]      verbosity     level of verbosity
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to associated linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_mumps_t *
cs_sles_mumps_create(int                          sym,
                     int                          verbosity,
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

  c->sym = sym;
  c->verbosity = verbosity;
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
    d = cs_sles_mumps_create(c->sym,
                             c->verbosity,
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

    sd->mumps->job = MUMPS_JOB_END;
    dmumps_c(sd->mumps);

    if (cs_glob_n_ranks == 1) {

      BFT_FREE(sd->mumps->irn);
      BFT_FREE(sd->mumps->jcn);
      BFT_FREE(sd->mumps->a);

    }
    else
      bft_error(__FILE__, __LINE__, 0, "%s: Not yet implemented", __func__);

    BFT_FREE(sd->mumps);
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

  cs_sles_mumps_t  *c = context;
  cs_sles_mumps_setup_t *sd = c->setup_data;

  if (sd == NULL) {
    BFT_MALLOC(c->setup_data, 1, cs_sles_mumps_setup_t);
    sd = c->setup_data;
  }

  int _verbosity = c->verbosity;
  if (_verbosity < 0)
    _verbosity = verbosity;

  BFT_MALLOC(sd->mumps, 1, DMUMPS_STRUC_C);

  /* Initialization step */

  sd->mumps->job = MUMPS_JOB_INIT;
  sd->mumps->par = 1;       /* all ranks are working */
  sd->mumps->sym = c->sym;

  /* Sanity checks */

  if (cs_matrix_get_diag_block_size(a)[0] > 1 ||
      cs_matrix_get_extra_diag_block_size(a)[0] > 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid matrix structure for MUMPS. No block requested.\n",
              __func__);

  const cs_halo_t  *halo = cs_matrix_get_halo(a);
  const cs_matrix_type_t  cs_mat_type = cs_matrix_get_type(a);
  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  if (halo != NULL)
    if (halo->n_transforms > 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: No periodicity requested with MUMPS.\n", __func__);

  /* Output */

  sd->mumps->ICNTL(1) = 6;      /* Error output: default value */

  if (_verbosity <= 0) {

    sd->mumps->ICNTL(2) = -1;   /* Rank statistics: default value */
    sd->mumps->ICNTL(3) = -1;   /* No global information printed */
    sd->mumps->ICNTL(4) = 0;    /* No message printed */
    sd->mumps->ICNTL(11) = 0;   /* No error analysis */

  }
  else {

    sd->mumps->ICNTL(2) = 0;    /* Rank statistics: default value */
    sd->mumps->ICNTL(3) = 6;    /* Global information: default value */
    sd->mumps->ICNTL(11) = 2;   /* Main error analysis */

    if (_verbosity == 1)
      sd->mumps->ICNTL(4) = 1;  /* Only error messages printed */
    else if (_verbosity == 2)
      sd->mumps->ICNTL(4) = 2;  /* Verbosity level: default value */
    else /* verbosity > 2 */
      sd->mumps->ICNTL(4) = 4;  /* All messages are printed */

  }

  if (cs_glob_n_ranks == 1) { /* sequential run */

    sd->mumps->ICNTL(5) = 0;    /* 0: assembled / 1: elemental */
    sd->mumps->ICNTL(18) = 0;   /* 0: centralized on rank 0 */
    sd->mumps->ICNTL(20) = 0;   /* 0: dense RHS on rank 0 */
    sd->mumps->ICNTL(21) = 0;   /* 0: dense solution array on rank 0 */

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

    sd->mumps->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(cs_glob_mpi_comm);
#else
    sd->mumps->comm_fortran = USE_COMM_WORLD; /* Not used in this case and set
                                                 to the default value given by
                                                 the MUMPS documentation */
#endif /* HAVE_MPI */

    /* First step: Initialization */

    dmumps_c(sd->mumps);

    /* Build the linear system for MUMPS */

    sd->mumps->n = (MUMPS_INT)n_rows;

    if (cs_mat_type == CS_MATRIX_MSR) {

      const cs_lnum_t  *a_row_idx, *a_col_ids;
      const cs_real_t  *d_val, *x_val;

      cs_matrix_get_msr_arrays(a, &a_row_idx, &a_col_ids, &d_val, &x_val);

      sd->mumps->nnz = (MUMPS_INT8)(n_rows + a_row_idx[n_rows]);

      BFT_MALLOC(sd->mumps->irn, sd->mumps->nnz, MUMPS_INT);
      BFT_MALLOC(sd->mumps->jcn, sd->mumps->nnz, MUMPS_INT);
      BFT_MALLOC(sd->mumps->a, sd->mumps->nnz, DMUMPS_REAL);

      /* Add diagonal entries */

      for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {

        sd->mumps->irn[row_id] = (MUMPS_INT)(row_id + 1);
        sd->mumps->jcn[row_id] = (MUMPS_INT)(row_id + 1);
        sd->mumps->a[row_id] = (DMUMPS_REAL)d_val[row_id];

      }

      /* Extra-diagonal entries */

      MUMPS_INT  *_irn = sd->mumps->irn + n_rows;
      MUMPS_INT  *_jcn = sd->mumps->jcn + n_rows;
      DMUMPS_REAL  *_a = sd->mumps->a + n_rows;

      if (c->sym == 1 || c->sym == 2) { /* Symmetric case */

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
      else { /* non-symmetric case */

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

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid matrix format", __func__);

  }
  else {

    bft_error(__FILE__, __LINE__, 0,
              " %s: Not implemented yet.\n", __func__);

  }

  if (c->setup_hook != NULL)
    c->setup_hook(c->hook_context, sd->mumps);

  /* Update return values */
  c->n_setups += 1;

  /* Second step: analysis and factorization */

  sd->mumps->job = MUMPS_JOB_ANALYSE_FACTO;
  dmumps_c(sd->mumps);

  if (cs_glob_rank_id < 1) {

    if (sd->mumps->INFOG(1) < 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n MUMPS feedback error code: INFOG(1)=%d, INFOG(2)=%d\n",
                    sd->mumps->INFOG(1), sd->mumps->INFOG(2));
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected during the anaylis/factorization step",
                __func__);
    }
    else {
      if (verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n MUMPS feedback code: INFOG(1)=%d, INFOG(2)=%d\n",
                      sd->mumps->INFOG(1), sd->mumps->INFOG(2));
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

  cs_sles_mumps_t  *c = context;
  cs_sles_mumps_setup_t  *sd = c->setup_data;

  if (sd == NULL) {
    cs_sles_mumps_setup(c, name, a, verbosity);
    sd = c->setup_data;
  }

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);

  /* Build the RHS */

  if (cs_glob_n_ranks == 1) { /* sequential run */

    /* Sanity checks */
    assert(n_rows == sd->mumps->n);
    assert(sizeof(cs_real_t) == sizeof(DMUMPS_REAL));

    sd->mumps->nrhs = 1;
    memcpy(vx, rhs, n_rows*sizeof(cs_real_t));
    sd->mumps->rhs = vx;

  }
  else {

    assert(cs_glob_n_ranks > 1);

    sd->mumps->nloc_rhs = sd->mumps->n;
    sd->mumps->lrhs_loc = sd->mumps->n;
    BFT_MALLOC(sd->mumps->rhs_loc, sd->mumps->n, DMUMPS_REAL);
    BFT_MALLOC(sd->mumps->irhs_loc, sd->mumps->n, MUMPS_INT);

  }

  /* Resolution */

  cs_fp_exception_disable_trap();

  sd->mumps->job = MUMPS_JOB_SOLVE;
  dmumps_c(sd->mumps);

  cs_fp_exception_restore_trap();

  /* Output */

  cs_sles_convergence_state_t cvg = CS_SLES_CONVERGED;
  if (sd->mumps->INFOG(1) < 0)
    cvg = CS_SLES_BREAKDOWN;

  *n_iter = 1;

  /* Scaled residual is stored in RINFOG(11) */
  *residue = sd->mumps->RINFOG(11);

  c->n_solves += 1;

  /* Solution is stored inside mumps->rhs */

  if (cs_glob_n_ranks == 1) {

    sd->mumps->rhs = NULL;

  }
  else {

    BFT_FREE(sd->mumps->rhs_loc);
    BFT_FREE(sd->mumps->irhs_loc);

    sd->mumps->rhs_loc = NULL;
    sd->mumps->irhs_loc = NULL;

  }

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

  switch(c->sym) {
  case 0:
    strncpy(sym_type_name, "unsymmetric", 31);
    break;
  case 1:
    strncpy(sym_type_name, "spd", 31);
    break;
  case 2:
    strncpy(sym_type_name, "symmetric", 31);
    break;
  default:
    snprintf(sym_type_name, 31, "%d", c->sym);
  }
  sym_type_name[31] = '\0';

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  "  Solver type:                       MUMPS\n"
                  "    Symm type:                         %s\n",
                  sym_type_name);

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
