/*============================================================================
 * Sparse Linear Equation Solver driver
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
#include "cs_blas.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_halo.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_util.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles.c

  \brief Sparse linear equation solver driver.

  The Sparse Linear Equation Solver subsystem is designed to allow
  both simple usage of built-in solvers, and plugging of solvers
  from external libraries.

  As the options associated with different solvers may be very varied,
  this subsystem is based on the use of a series of callback functions
  which may be associated with a given linear system (defined either
  by field id for quick and recurrent access, or by unique system
  name). It is possible to provide a default function so calls for
  the resolution of systems not specified before can be handled.

  To summarize, the functions here define a linear system solver
  driver, with the real work being done by functions bound to this model.
  The main intent is to help manage passing varied user options to the
  solvers, and handling consistency of logging.

  \enum cs_sles_convergence_state_t
  \brief Convergence status indicator.

  \var CS_SLES_DIVERGED
       The solver has diverged
  \var CS_SLES_BREAKDOWN
       The solver has broken down, and cannot make any more progress
  \var CS_SLES_MAX_ITERATION
       Maximum number of iterations has been reached, without reaching
       the required convergence
  \var CS_SLES_ITERATING
       The solver is iterating
  \var CS_SLES_CONVERGED
       The solver has converged

  \typedef  cs_sles_setup_t

  \brief  Function pointer for pre-resolution setup of a linear system's
          context.

  This setup may include building a multigrid hierarchy, or a preconditioner.

  Use of this type of function is optional: the context is expected to
  maintain state, so that if a cs_sles_solve_t function is called before a
  \ref cs_sles_setup_t function, the latter will be called automatically.

  \param[in, out]  context    pointer to solver context
  \param[in]       name       pointer to name of linear system
  \param[in]       a          matrix
  \param[in]       verbosity  associated verbosity

  \typedef  cs_sles_solve_t

  \brief  Function pointer for resolution of a linear system.

  If the associated cs_sles_setup_t function has not been called before
  this function, it will be called automatically.

  The solution context setup by this call (or that of the matching
  \ref cs_sles_setup_t function) will be maintained until the matching
  \ref cs_sles_free_t function is called.

  The matrix is not expected to change between successive calls, although
  the right hand side may. If the matrix changes, the associated
  \ref cs_sles_setup_t or \ref cs_sles_free_t function must be called between
  solves.

  The system is considered to have converged when
  residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.

  \param[in, out]  context        pointer to solver context
  \param[in]       name           pointer to name of linear system
  \param[in]       a              matrix
  \param[in]       verbosity      associated verbosity
  \param[in]       precision      solver precision
  \param[in]       r_norm         residual normalization
  \param[out]      n_iter         number of "equivalent" iterations
  \param[out]      residual       residual
  \param[in]       rhs            right hand side
  \param[out]      vx             system solution
  \param[in]       aux_size       number of elements in aux_vectors
  \param           aux_vectors    optional working area
                                  (internal allocation if NULL)

  \return  convergence status

  \typedef  cs_sles_free_t

  \brief  Function pointer for freeing of a linear system's context data.

  Note that this function should free resolution-related data, such as
  multigrid hierarchy, preconditioning, and any other temporary arrays or
  objects required for resolution, but should not free the whole context,
  as info used for logging (especially performance data) should be
  maintained.

  \param[in, out]  context  pointer to solver context.

  \typedef  cs_sles_log_t

  \brief  Function pointer for logging of linear solver history
          and performance data.

  This function will be called for each solver when \ref cs_sles_finalize
  is called.

  \param[in]  context   pointer to solver context
  \param[in]  log_type  log type

  \typedef  cs_sles_copy_t

  \brief  Function pointer for creation of a solver context based on
          the copy of another.

  The new context copies the settings of the copied context, but not
  its setup data and logged info, such as performance data.

  This type of function is optional, but enables associating different
  solvers to related systems (to differentiate logging) while using
  the same settings by default.

  \param[in]  context  context to copy

  \return  pointer to newly created context

  \typedef  cs_sles_destroy_t

  Function pointer for destruction of a linear system solver context.

  This function should free all context data, and will be called for each
  system when \ref cs_sles_finalize is called.

  \param[in, out]  context  pointer to solver context

  \typedef  cs_sles_error_handler_t

  \brief  Function pointer for handling of non-convergence when solving
          a linear system.

  Such a function is optional, and may be used for a variety of purposes,
  such as logging, postprocessing, re-trying with different parameters,
  aborting the run, or any combination thereof.

  An error handler may be  associated with a given solver context using
  \ref cs_sles_set_error_handler, in which case it will be called whenever
  convergence fails.

  \param[in, out]  sles           pointer to solver object
  \param[in]       status         convergence status
  \param[in]       a              matrix
  \param[in]       rhs            Right hand side
  \param[out]      vx             System solution

  \return  true if solve should be re-executed, false otherwise

  \typedef  cs_sles_define_t

  \brief  Function pointer for the default definition of a sparse
          linear equation solver

  The function may be associated using \ref cs_sles_set_default_define,
  so that it may provide a definition that will be used when
  \ref cs_sles_setup or \ref cs_sles_solve is used for a system for which
  no matching call to \ref cs_sles_define has been done.

  The function should call \ref cs_sles_define with arguments \p f_id
  and \p name, and appropriately chosen function pointers.

  A pointer to the matrix of the system to be solved is also provided,
  so that the corresponding information may be used to better choose
  defaults.

  \param[in]  f_id  associated field id, or < 0
  \param[in]  name  associated name if f_id < 0, or NULL
  \param[in]  a     matrix

  \typedef  cs_sles_verbosity_t

  \brief  Function pointer for the default definition of a sparse
          linear equation solver's verbosity

  The function may be associated using \ref cs_sles_set_default_verbosity, so
  that it may provide a definition that will be used when
  \ref cs_sles_find_or_add is called.

  \param[in]  f_id  associated field id, or < 0
  \param[in]  name  associated name if f_id < 0, or NULL

  \return  default verbosity value
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Postprocessing of linear system residuals */
/*-------------------------------------------*/

typedef struct {

  int                       writer_id;     /* writer id in case of
                                              postprocessing (0 for none) */

  cs_lnum_t                 n_rows;        /* number of rows */
  cs_lnum_t                 block_size;    /* size of block */

  cs_real_t                *row_residual;  /* residual */

} cs_sles_post_t;

/* Basic per linear system options and logging */
/*---------------------------------------------*/

struct _cs_sles_t {

  int                       n_calls;       /* Number of setup or solve
                                              calls */

  int                       n_no_op;       /* Number of solves with immediate
                                              exit */
  bool                      allow_no_op;   /* Allow immediate exit if RHS
                                              small relative to residual norm */

  int                       f_id;          /* matching field id, or < 0 */

  const char               *name;          /* name if f_id < 0, or NULL */
  char                     *_name;         /* private name if f_id < 0,
                                              or NULL */

  int                       verbosity;     /* verbosity level */

  int                       type_id;       /* id of solver type */
  void                     *context;       /* solver context
                                              (options, state, logging) */

  cs_sles_setup_t          *setup_func;    /* solver setup function */
  cs_sles_solve_t          *solve_func;    /* solve function */
  cs_sles_free_t           *free_func;     /* free setup function */

  cs_sles_log_t            *log_func;      /* logging function */

  cs_sles_copy_t           *copy_func;     /* copy function */

  cs_sles_destroy_t        *destroy_func;  /* destruction function */

  cs_sles_error_handler_t  *error_func;    /* error handler */

  cs_sles_post_t           *post_info;     /* postprocessing info */

};

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Type name map */

static cs_map_name_to_id_t  *_type_name_map = NULL;

/* Pointers to default definitions */

static cs_sles_define_t *_cs_sles_define_default = NULL;
static cs_sles_verbosity_t *_cs_sles_default_verbosity = NULL;

/* Current and maximum number of systems respectively defined by field id,
   by name, or redefined after use */

static int _cs_sles_n_systems[3] = {0, 0, 0};
static int _cs_sles_n_max_systems[3] = {0, 0, 0};

/* Arrays of settings and status for linear systems */
static cs_sles_t **_cs_sles_systems[3] = {NULL, NULL, NULL};

/* Timer statistics */

static cs_timer_counter_t   _sles_t_tot;     /* Total time in linear solvers */
static int _sles_stat_id = -1;

/* Threshold for the detection of an immediate exit */

static double _cs_sles_epzero = 1e-12;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create new sparse linear equation solver structure.
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   pointer to associated linear system object.
 *----------------------------------------------------------------------------*/

static cs_sles_t *
_sles_create(int          f_id,
             const char  *name)
{
  cs_sles_t *sles;

  BFT_MALLOC(sles, 1, cs_sles_t);

  sles->f_id = f_id;

  if (f_id < 0 && name != NULL) {
    BFT_MALLOC(sles->_name, strlen(name) + 1, char);
    strcpy(sles->_name, name);
  }
  else
    sles->_name = NULL;

  if (_cs_sles_default_verbosity != NULL)
    sles->verbosity = _cs_sles_default_verbosity(f_id, name);
  else
    sles->verbosity = 0;

  if (_type_name_map == NULL)
    _type_name_map = cs_map_name_to_id_create();
  sles->type_id = cs_map_name_to_id(_type_name_map, "<undefined>");

  sles->name = sles->_name;

  sles->context = NULL;
  sles->setup_func = NULL;
  sles->solve_func = NULL;
  sles->free_func = NULL;
  sles->log_func = NULL;
  sles->copy_func = NULL;
  sles->destroy_func = NULL;
  sles->error_func = NULL;

  sles->n_calls = 0;
  sles->n_no_op = 0;
  sles->allow_no_op = false;

  sles->post_info = NULL;

  return sles;
}

/*----------------------------------------------------------------------------
 * Return pointer to linear system object, based on matching field id.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   f_id <-- associated field id
 *
 * returns:
 *   pointer to linear system object
 *----------------------------------------------------------------------------*/

static cs_sles_t *
_find_or_add_system_by_f_id(int  f_id)
{
  assert(f_id >= 0);

  /* Return system id already defined */
  if (f_id < _cs_sles_n_max_systems[0]) {
    if (_cs_sles_systems[0][f_id] != NULL)
      return _cs_sles_systems[0][f_id];
  }

  /* Increase size of systems array if required */
  else {
    int i = _cs_sles_n_max_systems[0];

    if (_cs_sles_n_max_systems[0] == 0)
      _cs_sles_n_max_systems[0] = 1;
    while (_cs_sles_n_max_systems[0] <= f_id)
      _cs_sles_n_max_systems[0] *= 2;
    BFT_REALLOC(_cs_sles_systems[0],
                _cs_sles_n_max_systems[0],
                cs_sles_t *);

    for (int j = i; j < _cs_sles_n_max_systems[0]; j++)
      _cs_sles_systems[0][j] = NULL;
  }

  /* If we have not returned yet, we need to add a new system to the array */

  _cs_sles_systems[0][f_id] = _sles_create(f_id, NULL);
  _cs_sles_n_systems[0] += 1;

  return _cs_sles_systems[0][f_id];
}

/*----------------------------------------------------------------------------
 * Return pointer to linear system object, based on system name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   name <-- system name
 *
 * returns:
 *   pointer to linear system object
 *----------------------------------------------------------------------------*/

static cs_sles_t *
_find_or_add_system_by_name(const char  *name)
{
  int ii, start_id, end_id, mid_id;
  int cmp_ret = 1;

  /* Use binary search to find system */

  start_id = 0;
  end_id = _cs_sles_n_systems[1] - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp((_cs_sles_systems[1][mid_id])->name, name);
    if (cmp_ret < 0)
      start_id = mid_id + 1;
    else if (cmp_ret > 0)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  /* If found, return */

  if (cmp_ret == 0)
    return _cs_sles_systems[1][mid_id];

  /* Reallocate global array if necessary */

  if (_cs_sles_n_systems[1] >= _cs_sles_n_max_systems[1]) {
    int i = _cs_sles_n_max_systems[1];

    if (_cs_sles_n_max_systems[1] == 0)
      _cs_sles_n_max_systems[1] = 1;
    _cs_sles_n_max_systems[1] *= 2;
    BFT_REALLOC(_cs_sles_systems[1],
                _cs_sles_n_max_systems[1],
                cs_sles_t*);

    for (int j = i; j < _cs_sles_n_max_systems[1]; j++)
      _cs_sles_systems[1][j] = NULL;
  }

  /* Insert in sorted list */

  for (ii = _cs_sles_n_systems[1]; ii > mid_id; ii--)
    _cs_sles_systems[1][ii] = _cs_sles_systems[1][ii - 1];

  _cs_sles_systems[1][mid_id] = _sles_create(-1, name);
  _cs_sles_n_systems[1] += 1;

  return _cs_sles_systems[1][mid_id];
}

/*----------------------------------------------------------------------------
 * Copy linear system info whose definition has changed.
 *
 * This is intended to maintain logging info, and is useful only
 * if the solver has been called using the older definition.
 *
 * parameters:
 *   s <-> pointer to linear system info
 *----------------------------------------------------------------------------*/

static void
_save_system_info(cs_sles_t  *s)
{
  assert(s != NULL);

  /* Resize array if needed */

  int i = _cs_sles_n_systems[2];

  if (i >= _cs_sles_n_max_systems[2]) {

    if (_cs_sles_n_max_systems[2] == 0)
      _cs_sles_n_max_systems[2] = 1;
    _cs_sles_n_max_systems[2] *=2;
    BFT_REALLOC(_cs_sles_systems[2],
                _cs_sles_n_max_systems[2],
                cs_sles_t *);

    for (int j = i; j < _cs_sles_n_max_systems[2]; j++)
      _cs_sles_systems[2][j] = NULL;

  }

  /* Ensure no extra data is maintained for old system */

  if (s->free_func != NULL)
    s->free_func(s->context);

  /* Save other context and options */

  cs_sles_t *s_old;
  BFT_MALLOC(s_old, 1, cs_sles_t);
  memcpy(s_old, s, sizeof(cs_sles_t));

  s_old->_name = NULL; /* still points to new name */
  s->context = NULL;   /* old context now only available through s_old */

  _cs_sles_systems[2][i] = s_old;

  _cs_sles_n_systems[2] += 1;
}

/*----------------------------------------------------------------------------
 * Compute per-cell residual for Ax = b.
 *
 * parameters:
 *   n_vals           <-- Number of values
 *   a                <-- Linear equation matrix
 *   rhs              <-- Right hand side
 *   vx               <-> Current system solution
 *   res              --> Residual
 *----------------------------------------------------------------------------*/

static void
_residual(cs_lnum_t            n_vals,
          const cs_matrix_t   *a,
          const cs_real_t      rhs[],
          cs_real_t            vx[],
          cs_real_t            res[])
{
#if defined(HAVE_ACCEL)

  bool ddp_vx = false, ddp_res = false;

  if (cs_matrix_get_alloc_mode(a) > CS_ALLOC_HOST) {
    cs_lnum_t _n_vals =   cs_matrix_get_n_columns(a)
                        * cs_matrix_get_diag_block_size(a);
    if (cs_check_device_ptr(vx) == CS_ALLOC_HOST) {
      cs_associate_device_ptr(vx, _n_vals, sizeof(cs_real_t));
      ddp_vx = true;
    }
    if (cs_check_device_ptr(res) == CS_ALLOC_HOST) {
      cs_associate_device_ptr(res, _n_vals, sizeof(cs_real_t));
      ddp_res = true;
    }
  }

#endif

  cs_matrix_vector_multiply(a, vx, res);

# pragma omp parallel for if(n_vals > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_vals; ii++)
    res[ii] = fabs(res[ii] - rhs[ii]);

#if defined(HAVE_ACCEL)

  if (ddp_vx)
    cs_disassociate_device_ptr(vx);
  if (ddp_res)
    cs_disassociate_device_ptr(res);

#endif
}

/*----------------------------------------------------------------------------
 * Test if a general sparse linear system needs solving or if the right-hand
 * side is already zero within convergence criteria.
 *
 * The computed residual is also updated;
 *
 * parameters:
 *   name      <-- name of the associated system
 *   a         <-- pointer to matrix
 *   verbosity <-- verbosity level
 *   precision <-- solver precision
 *   r_norm    <-- residual normalization
 *   residual  <-> residual
 *   vx        <-- initial solution
 *   rhs       <-- right hand side
 *
 * returns:
 *   1 if solving is required, 0 if the rhs is already zero within tolerance
 *   criteria (precision of residual normalization)
 *----------------------------------------------------------------------------*/

static int
_needs_solving(const  char        *name,
               const cs_matrix_t  *a,
               int                 verbosity,
               double              precision,
               double              r_norm,
               double             *residual,
               const cs_real_t    *vx,
               const cs_real_t    *rhs)
{
  int retval = 1;

  /* Initialize residual, check for immediate return */

  const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a) * diag_block_size;

  double r[2] = {
    cs_dot_xx(n_rows, rhs),
    cs_dot_xx(n_rows, vx)
  };
  cs_parall_sum(2, CS_DOUBLE, r);

  /* If the initial solution is "true" zero (increment mode), we can determine
     convergence without resorting to a matrix-vector product */

  if (r[1] < 1e-60) {

    double _precision = CS_MIN(_cs_sles_epzero, /* prefer to err on the side */
                               precision);      /* of caution... */

    *residual = sqrt(r[0]);

    if (r_norm <= _cs_sles_epzero)
      retval = 0;
    else if (*residual/r_norm <= _precision)
      retval = 0;

    if (retval == 0 && verbosity > 1)
      bft_printf(_("[%s]:\n"
                   "  immediate exit; r_norm = %11.4e, residual = %11.4e\n"),
                 name, r_norm, *residual);

  }
  else
    *residual = HUGE_VAL; /* actually unknown, since we did not multiply
                            by A (we might as well enter the solver,
                            and expect to have vx = 0 most of the time) */

  return retval;
}

/*----------------------------------------------------------------------------
 * Output post-processing data for failed system convergence.
 *
 * parameters:
 *   n_vals        <-- Size of val and val_type array
 *   val           <-> Values to post-process (set to 0 on output if not
 *                     normal floating-point values)
 *   val_type      --> 0: normal values, 1: infinite, 2: Nan
 *
 * returns:
 *   number of non-normal values
 *----------------------------------------------------------------------------*/

static size_t
_value_type(size_t     n_vals,
            cs_real_t  val[],
            cs_real_t  val_type[])
{
  size_t ii;
  size_t retval = 0;

#if (__STDC_VERSION__ >= 199901L)

  for (ii = 0; ii < n_vals; ii++) {

    int v_type = fpclassify(val[ii]);

    if (v_type == FP_INFINITE) {
      val[ii] = 0.;
      val_type[ii] = 1;
      retval += 1;
    }

    else if (v_type == FP_NAN) {
      val[ii] = 0.;
      val_type[ii] = 2;
      retval += 1;
    }

    else if (val[ii] > 1.e38 || val[ii] < -1.e38) {
      val[ii] = 0.;
      val_type[ii] = 1;
      retval += 1;
    }

    else
      val_type[ii] = 0;
  }

#else

  for (ii = 0; ii < n_vals; ii++) {

    if (val[ii] != val[ii]) { /* Test for NaN with IEEE 754 arithmetic */
      val[ii] = 0.;
      val_type[ii] = 2;
      retval += 1;
    }

    else if (val[ii] > 1.e38 || val[ii] < -1.e38) {
      val[ii] = 0.;
      val_type[ii] = 1;
      retval += 1;
    }

    else
      val_type[ii] = 0;
  }

#endif

  return retval;
}

/*----------------------------------------------------------------------------
 * Ensure array for postprocessing output of residuals is allocated
 *
 * parameters
 *   sles <-> pointer to solver object
 *   a    <-- matrix
 *----------------------------------------------------------------------------*/

static void
_ensure_alloc_post(cs_sles_t          *sles,
                   const cs_matrix_t  *a)
{
  if (sles->post_info != NULL) {
    const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);
    const cs_lnum_t n_vals = cs_matrix_get_n_columns(a) * diag_block_size;

    sles->post_info->n_rows = cs_matrix_get_n_rows(a);
    sles->post_info->block_size = diag_block_size;

    BFT_REALLOC(sles->post_info->row_residual, n_vals, cs_real_t);
  }
}

/*----------------------------------------------------------------------------
 * Post process the residual for a given linear equation solver
 *
 * parameters:
 *   sles_p <-- void pointer to sparse linear equation solver context
 *   ts     <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_post_function(void                  *sles_p,
               const cs_time_step_t  *ts)
{
  CS_UNUSED(ts);

  cs_sles_t *sles = sles_p;

  cs_sles_post_t *sp = sles->post_info;

  assert(sp != NULL);

  const cs_mesh_t *mesh = cs_glob_mesh;

  int mesh_id = CS_POST_MESH_VOLUME;

  char base_name[32], val_name[32];

  /* Check for mesh location */

  int location_id = 0, flag = 0;
  if (sp->n_rows == mesh->n_cells)
    location_id = CS_MESH_LOCATION_CELLS;
  else if (sp->n_rows == mesh->n_vertices)
    location_id = CS_MESH_LOCATION_VERTICES;

  flag = location_id;
  cs_parall_max(1, CS_INT_TYPE, &flag);
  if (flag == location_id)
    flag = 0;
  else
    flag =1;
  cs_parall_max(1, CS_INT_TYPE, &flag);
  if (flag != 0)
    return;

  strcpy(base_name, "Residual");

  const char *name = cs_sles_get_name(sles);

  if (strlen(name) + strlen(base_name) < 31) {
    strcpy(val_name, base_name);
    strcat(val_name, "_");
    strcat(val_name, name);
  }
  else {
    strncpy(val_name, base_name, 31);
    val_name[31] = '\0';
  }

  cs_sles_post_output_var(val_name,
                          mesh_id,
                          location_id,
                          sp->writer_id,
                          sp->block_size,
                          sp->row_residual);

  BFT_FREE(sp->row_residual);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the threshold value used in the detection of immediate exit
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_epzero(double  new_value)
{
  _cs_sles_epzero = new_value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the current threshold value used in the detection of immediate
 *        exit
 *
 * \return the value of the threshold
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_get_epzero(void)
{
  return _cs_sles_epzero;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize sparse linear equation solver API.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_initialize(void)
{
  CS_TIMER_COUNTER_INIT(_sles_t_tot);

  int stats_root = cs_timer_stats_id_by_name("operations");

  if (stats_root > -1) {
    _sles_stat_id = cs_timer_stats_create("operations",
                                          "linear_solvers",
                                          "linear solvers");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize sparse linear equation solver API.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_finalize(void)
{
  for (int i = 0; i < 3; i++) {

    for (int j = 0; j < _cs_sles_n_max_systems[i]; j++) {

      if (_cs_sles_systems[i][j] != NULL) {
        cs_sles_t *sles = _cs_sles_systems[i][j];
        if (sles->free_func != NULL)
          sles->free_func(sles->context);
        if (sles->destroy_func != NULL)
          sles->destroy_func(&(sles->context));
        if (sles->post_info != NULL) {
          BFT_FREE(sles->post_info->row_residual);
          BFT_FREE(sles->post_info);
        }
        BFT_FREE(sles->_name);
        BFT_FREE(_cs_sles_systems[i][j]);
      }

    }

    BFT_FREE(_cs_sles_systems[i]);
    _cs_sles_n_max_systems[i] = 0;
    _cs_sles_n_systems[i] = 0;

  }

  cs_map_name_to_id_destroy(&_type_name_map);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info
 *
 * \param[in]  log_type  log type (CS_LOG_SETUP or CS_LOG_PERFORMANCE)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_log(cs_log_t  log_type)
{
  int n_tot_systems = 0;
  for (int i = 0; i < 3; i++)
    n_tot_systems += _cs_sles_n_systems[i];

  if (n_tot_systems < 1)
    return;

  int log_order[] = {2, 0, 1}; /* log previous setups first, then fields,
                                  then others */

  if (log_type == CS_LOG_PERFORMANCE)
    cs_log_printf
      (log_type,
       _("\n"
         "Total elapsed time for linear equation system solvers:  %.3f s\n"),
       _sles_t_tot.nsec*1e-9);

  else if (log_type == CS_LOG_SETUP) {

    const char  header[] = "Linear solver options for all systems";
    size_t l = cs_log_strlen(header);
    char ul[128];

    l = CS_MIN(l, 127);
    for (size_t ll = 0; ll < l; ll++)
      ul[ll] = '-';
    ul[l] = '\0';

    cs_log_printf(log_type, "\n%s\n", header);
    cs_log_printf(log_type, "%s\n\n", ul);
    cs_log_printf(log_type,
                  "Immediate exit threshold value: %5.2e\n", _cs_sles_epzero);

  }

  const char *option_category[]
    = {N_("Linear solver options modified during run (previous values)"),
       N_("Linear solver options for fields"),
       N_("Linear solver options for other systems")};

  const char *perf_category[]
    = {N_("Linear solver performance with previous options"),
       N_("Linear solver performance for fields"),
       N_("Linear solver performance for other systems")};

  for (int i = 0; i < 3; i++) {

    int j = log_order[i];

    /* Print heading (underlined) line */

    if (_cs_sles_n_systems[j] > 0) {

      size_t l = 0;
      switch(log_type) {
      case CS_LOG_SETUP:
        l = cs_log_strlen(_(option_category[i]));
        cs_log_printf(log_type, "\n%s\n", _(option_category[i]));
        break;
      case CS_LOG_PERFORMANCE:
        l = cs_log_strlen(_(perf_category[i]));
        cs_log_printf(log_type, "\n%s\n", _(perf_category[i]));
        break;
      default:
        break;
      }

      char ul[128];
      l = CS_MIN(l, 127);
      for (size_t ll = 0; ll < l; ll++)
        ul[ll] = '-';
      ul[l] = '\0';
      cs_log_printf(log_type, "%s\n", ul);

    }

    /* Logging for each system */

    for (int k = 0; k < _cs_sles_n_max_systems[j]; k++) {
      cs_sles_t *sles = _cs_sles_systems[j][k];

      if (sles == NULL) continue;

      if (sles->log_func != NULL) {

        const char *name = cs_sles_base_name(sles->f_id, sles->name);

        switch(log_type) {

        case CS_LOG_SETUP:
          if (sles->f_id > -1)
            cs_log_printf
              (log_type,
               _("\n"
                 "Linear solver options for \"%s\" (field id %d)\n"),
               name, sles->f_id);
          else
            cs_log_printf
              (log_type,
               _("\n"
                 "Linear solver options for \"%s\"\n"),
               name);
          break;

        case CS_LOG_PERFORMANCE:
          if (sles->f_id > -1)
            cs_log_printf
              (log_type,
               _("\n"
                 "Summary of resolutions for \"%s\" (field id %d)\n"),
               name, sles->f_id);
          else
            cs_log_printf
              (log_type,
               _("\n"
                 "Summary of resolutions for \"%s\"\n"),
               name);
          break;

        default:
          break;
        }

        sles->log_func(sles->context, log_type);

        switch(log_type) {

        case CS_LOG_SETUP:
          cs_log_printf
            (log_type,
             _("  Verbosity: %d\n"), sles->verbosity);
          if (sles->post_info != NULL)
            cs_log_printf
              (log_type,
               _("  Residual postprocessing writer id: %d\n"),
               sles->post_info->writer_id);
          break;

        case CS_LOG_PERFORMANCE:
          if (sles->n_no_op > 0)
            cs_log_printf
              (log_type,
               _("\n"
                 "  Number of immediate solve exits: %d\n"), sles->n_no_op);
          break;

        default:
          break;
        }

      }
    }

  }

  cs_log_printf(log_type, "\n");
  cs_log_separator(log_type);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to linear system object, based on matching field id or
 *        system name.
 *
 * If this system did not previously exist, NULL is returned.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 *
 * \return  pointer to associated linear system object, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_find(int          f_id,
             const char  *name)
{
  cs_sles_t *retval = NULL;

  if (f_id >= 0) {
    if (f_id < _cs_sles_n_max_systems[0]) {
      if (_cs_sles_systems[0][f_id] != NULL) {
        retval = _cs_sles_systems[0][f_id];
        /* Check for masked ("pushed") definition */
        if (retval->name != NULL)
          retval = cs_sles_find(-1, retval->name);
      }
    }
  }

  else if (name != NULL) {

    int start_id, end_id, mid_id;
    int cmp_ret = 1;

    /* Use binary search to find system */

    start_id = 0;
    end_id = _cs_sles_n_systems[1] - 1;
    mid_id = start_id + ((end_id -start_id) / 2);

    while (start_id <= end_id) {
      cmp_ret = strcmp((_cs_sles_systems[1][mid_id])->name, name);
      if (cmp_ret < 0)
        start_id = mid_id + 1;
      else if (cmp_ret > 0)
        end_id = mid_id - 1;
      else
        break;
      mid_id = start_id + ((end_id -start_id) / 2);
    }

    if (cmp_ret == 0)
      retval = _cs_sles_systems[1][mid_id];
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to linear system object, based on matching field id or
 *        system name.
 *
 * If this system did not previously exist, it is created and added to
 * to the list of "known" systems. In this case, it will be usable
 * only if cs_sles_define() is called for the same field id and name
 * (in which case calling the present function is redundant), or if
 * cs_sles_set_sefault_define() has been previously used to define
 * the default solver policy.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 *
 * \return  pointer to associated linear system object, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_find_or_add(int          f_id,
                    const char  *name)
{
  cs_sles_t *retval = NULL;

  if (f_id >= 0) {
    retval = _find_or_add_system_by_f_id(f_id);
    /* Check for masked ("pushed") definition */
    if (retval->name != NULL)
      retval = _find_or_add_system_by_name(retval->name);
  }
  else
    retval = _find_or_add_system_by_name(name);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Temporarily replace field id with name for matching calls
 * to \ref cs_sles_setup, \ref cs_sles_solve, \ref cs_sles_free, and other
 * operations involving access through a field id.
 *
 * This function is provided to allow some peculiar calling sequences,
 * in which \ref cs_equation_iterative_solve_scalar is called with a given
 * field id, but specific solver options must still be set.
 * In the future, a cleaner method to handle those exceptional cases
 * would be preferred. As such, only a stack depth of 1 is allowed.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_push(int          f_id,
             const char  *name)
{
  if (f_id < 0)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s must be called only for an actual field, with id >=0, not %d.",
       __func__, f_id);

  cs_sles_t *retval = cs_sles_find_or_add(f_id, NULL);

  if (retval->name != NULL)
    bft_error
      (__FILE__, __LINE__, 0,
       _("cs_sles_push() only allows a stack of depth 1:\n"
         "  it  may not be called multiple times for a given field (id %d)\n"
         "  without calling cs_sles_pop between those calls."), f_id);
  else {
    assert(retval->_name == NULL);
    BFT_MALLOC(retval->_name, strlen(name) + 1, char);
    strcpy(retval->_name, name);
    retval->name = retval->_name;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Restore behavior temporarily modified by \ref cs_sles_push.
 *
 * \deprecated This function matches \ref cs_sles_push, which is deprecated.
 *
 * \param[in]  f_id  associated field id, or < 0
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pop(int  f_id)
{
  if (f_id < 0)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s must be called only for an actual field, with id >=0, not %d.",
       __func__, f_id);

  cs_sles_t *retval = _find_or_add_system_by_f_id(f_id);

  retval->name = NULL;
  BFT_FREE(retval->_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define sparse linear equation solver for a given field or
 *        equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * The context pointer is used to point to a structure adapted to the function
 * pointers given here, and combined with those functions, allows using
 * both built-in, external, or user-defined solvers.
 *
 * It is recommended the context type name provided here directly relate
 * to the associated structure type (for example, "cs_sles_it_t" or
 * "cs_multigrid_t").
 *
 * \param[in]       f_id          associated field id, or < 0
 * \param[in]       name          associated name if f_id < 0, or NULL
 * \param[in, out]  context       pointer to solver context management
 *                                structure (cs_sles subsystem becomes owner)
 * \param[in]       type_name     context structure or object type name
 * \param[in]       setup_func    pointer to system setup function
 * \param[in]       solve_func    pointer to system solution function (also
 *                                calls setup_func if not done yet or free_func
 *                                called since last solve)
 * \param[in]       free_func     pointer function freeing system setup
 * \param[in]       log_func      pointer to system info logging function
                                  (optional, but recommended)
 * \param[in]       copy_func     pointer to settings copy function (optional)
 * \param[in]       destroy_func  pointer to function destroying solver context
 *                                (called with \ref cs_sles_finalize or with a
 *                                new call to this function for the same system)
 *
 * \return  pointer to associated linear system object
 */
/*----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_define(int                 f_id,
               const char         *name,
               void               *context,
               const char         *type_name,
               cs_sles_setup_t    *setup_func,
               cs_sles_solve_t    *solve_func,
               cs_sles_free_t     *free_func,
               cs_sles_log_t      *log_func,
               cs_sles_copy_t     *copy_func,
               cs_sles_destroy_t  *destroy_func)
{
  cs_sles_t * sles = cs_sles_find_or_add(f_id, name);

  /* Check if system was previously defined and used,
     and save info for future logging in this case */

  if (sles->context != NULL) {
    if (sles->n_calls > 0  && sles->log_func != NULL)
      _save_system_info(sles);
    else if (sles->destroy_func != NULL)
      sles->destroy_func(&(sles->context));
  }

  if (type_name != NULL)
    sles->type_id = cs_map_name_to_id(_type_name_map, type_name);

  /* Now define options */

  sles->context = context;
  sles->setup_func = setup_func;
  sles->solve_func = solve_func;
  sles->free_func = free_func;
  sles->log_func = log_func;
  sles->copy_func = copy_func;
  sles->destroy_func = destroy_func;

  return sles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the verbosity for a given linear equation solver.
 *
 * This verbosity will be used by cs_sles_setup and cs_sles_solve.
 *
 * By default, the verbosity is set to 0, or the value returned by the
 * function set with cs_sles_set_default_define().
 *
 * \param[in, out]  sles       pointer to solver object
 * \param[in]       verbosity  verbosity level
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_verbosity(cs_sles_t  *sles,
                      int         verbosity)
{
  sles->verbosity = verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the verbosity for a given linear equation solver.
 *
 * This verbosity will be used by cs_sles_setup and cs_sles_solve.
 *
 * \param[in, out]  sles       pointer to solver object
 *
 * \return  verbosity level
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_verbosity(cs_sles_t  *sles)
{
  return sles->verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate postprocessing output for a given linear equation solver.
 *
 * This allows the output of the residual at the end of each solution
 * series, using a single postprocessing writer.
 * By default, no output is activated.
 *
 * \param[in, out]  sles       pointer to solver object
 * \param[in]       writer_id  id of the writer
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_post_output(cs_sles_t  *sles,
                        int         writer_id)
{
  if (sles->n_calls > 0)
    return;

  if (sles->post_info == NULL)
    cs_post_add_time_dep_output(_post_function, (void *)sles);

  BFT_REALLOC(sles->post_info, 1, cs_sles_post_t);
  sles->post_info->writer_id = writer_id;
  sles->post_info->n_rows = 0;
  sles->post_info->block_size = 0;
  sles->post_info->row_residual = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of the associated writer if postprocessing output
 *        is active for a given linear equation solver.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  id od associated writer, or 0
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_post_output(cs_sles_t  *sles)
{
  int retval = 0;

  if (sles->post_info != NULL)
    retval = sles->post_info->writer_id;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return type name of solver context.
 *
 * The returned string is intended to help determine which type is associated
 * with the void * pointer returned by \ref cs_sles_get_context for a given
 * solver definition, so as to be able to call additional specific functions
 * beyond the generic functions assigned to a cs_sles_t object.
 *
 * If no type_name string was associated to the solver upon its definition by
 * \ref cs_sles_define, or it has not been defined yet, the string returned
 * is "<undefined>". It is recommended the type name provided
 * \ref cs_sles_define directly relate to the associated structure type
 * (for example, "cs_sles_it_t" or "cs_multigrid_t").
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  pointer to linear system solver specific type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_get_type(cs_sles_t  *sles)
{
  return cs_map_name_to_id_reverse(_type_name_map, sles->type_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to solver context structure pointer.
 *
 * The context structure depends on the type of solver used, which may in
 * turn be determined by the string returned by cs_sles_get_type().
 * If may be used by appropriate functions specific to that type.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  pointer to solver-specific linear system info and context
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_get_context(cs_sles_t  *sles)
{
  return sles->context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return field id associated with a given sparse linear equation solver.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  associated field id (or -1 if defined by name)
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_f_id(const cs_sles_t  *sles)
{
  return sles->f_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return name associated with a given sparse linear equation solver.
 *
 * This is simply a utility function which will return its name argument
 * if f_id < 0, and the associated field's name or label otherwise.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  pointer to associated linear system object name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_get_name(const cs_sles_t  *sles)
{
  return cs_sles_base_name(sles->f_id, sles->name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query if immediate_return ("no-op") is allowed when initial
 * guess is zero (solve by increments) and the RHS is already zero within the
 * normalized tolerance criteria.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  true if immediate return is allowed, false if at least one
 *          iteration is required
 */
/*----------------------------------------------------------------------------*/

bool
cs_sles_get_allow_no_op(const cs_sles_t  *sles)
{
  return sles->allow_no_op;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if immediate_return ("no-op") is allowed when initial
 * guess is zero (solve by increments) and the RHS is already zero within the
 * normalized tolerance criteria.
 *
 * \param[in, out]  sles         pointer to solver object
 * \param[in]       allow_no_op  true if immediate return is allowed,
 *                               false if at least one iteration is required
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_allow_no_op(cs_sles_t  *sles,
                        bool        allow_no_op)
{
  sles->allow_no_op = allow_no_op;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup sparse linear equation solver.
 *
 * Use of this function is optional: if a \ref cs_sles_solve is called
 * for the same system before this function is called, the latter will be
 * called automatically.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * \param[in, out]  sles  pointer to solver object
 * \param[in]       a     matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup(cs_sles_t          *sles,
              const cs_matrix_t  *a)
{
  cs_timer_t t0 = cs_timer_time();

  if (sles->context == NULL)
    _cs_sles_define_default(sles->f_id, sles->name, a);

  int t_top_id = cs_timer_stats_switch(_sles_stat_id);

  sles->n_calls += 1;

  if (sles->setup_func != NULL) {
    const char  *sles_name = cs_sles_base_name(sles->f_id, sles->name);
    sles->setup_func(sles->context, sles_name, a, sles->verbosity);
  }

  /* Prepare residual postprocessing if required */

  if (sles->post_info != NULL) {
    _ensure_alloc_post(sles, a);
    const cs_lnum_t n_vals
      = cs_matrix_get_n_columns(a) * sles->post_info->block_size;
    cs_real_t *r = sles->post_info->row_residual;
#   pragma omp parallel for if(n_vals > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vals; i++)
      r[i] = 0;
  }

  cs_timer_stats_switch(t_top_id);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_sles_t_tot, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief General sparse linear system resolution.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * Note that if \ref cs_sles_setup was previously called for this
 * system, and \ref cs_sles_free has not been called since, the matrix
 * provided should be the same. The optional separation between the
 * two stages is intended to allow amortizing the cost of setup
 * over multiple solutions.
 *
 * The system is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 *
 * \param[in, out]  sles           pointer to solver object
 * \param[in]       a              matrix
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       size of aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if NULL)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve(cs_sles_t           *sles,
              const cs_matrix_t   *a,
              double               precision,
              double               r_norm,
              int                 *n_iter,
              double              *residual,
              const cs_real_t     *rhs,
              cs_real_t           *vx,
              size_t               aux_size,
              void                *aux_vectors)
{
  cs_timer_t t0 = cs_timer_time();

  if (sles->context == NULL)
    _cs_sles_define_default(sles->f_id, sles->name, a);

  int t_top_id = cs_timer_stats_switch(_sles_stat_id);

  sles->n_calls += 1;

  assert(sles->solve_func != NULL);

  const char  *sles_name = cs_sles_base_name(sles->f_id, sles->name);

#if 0
  /* Dump linear system to file (for experimenting with external tools) */
  cs_matrix_dump_linear_system(a, rhs, sles_name);
#endif

  cs_sles_convergence_state_t state;
  bool do_solve = true;

  /* Even if we normally require at least entering the linear equation solver,
     if the residual normalization is really zero, we probably have zero initial
     solution and RHS already, so check for that case, rather than enter
     solvers whose convergence test may fail in these conditions. */

  if (sles->allow_no_op || r_norm <= 0.) {
    do_solve = _needs_solving(sles_name,
                              a,
                              sles->verbosity,
                              precision,
                              r_norm,
                              residual,
                              vx,
                              rhs);

    if (! do_solve) {
      sles->n_no_op += 1;
      *n_iter = 0;
      state = CS_SLES_CONVERGED;
    }
  }

  while (do_solve) {

    state = sles->solve_func(sles->context,
                             sles_name,
                             a,
                             sles->verbosity,
                             precision,
                             r_norm,
                             n_iter,
                             residual,
                             rhs,
                             vx,
                             aux_size,
                             aux_vectors);

    if (state < CS_SLES_ITERATING && sles->error_func != NULL)
      do_solve = sles->error_func(sles,
                                  state,
                                  a,
                                  rhs,
                                  vx);
    else
      do_solve = false;

  }

  /* Prepare postprocessing if needed */

  if (sles->post_info != NULL) {
    _ensure_alloc_post(sles, a);
    const cs_lnum_t n_vals
      = sles->post_info->n_rows * sles->post_info->block_size;
    _residual(n_vals,
              a,
              rhs,
              vx,
              sles->post_info->row_residual);
  }

  /* Check error */

  if (sles->verbosity > 1) {
    const cs_lnum_t block_size = cs_matrix_get_diag_block_size(a);
    const cs_lnum_t n_vals_ext = cs_matrix_get_n_columns(a) * block_size;
    const cs_lnum_t n_vals = cs_matrix_get_n_rows(a) * block_size;

    cs_real_t *resr = NULL;
    BFT_MALLOC(resr, n_vals_ext, cs_real_t);

    _residual(n_vals, a, rhs, vx, resr);
    cs_real_t rsd = sqrt(cs_gdot(n_vals, resr, resr));

    bft_printf
      ("# residual[%s] = %g (%g * required, precision %g, normalization %g)\n",
       sles_name, rsd, rsd/(precision*r_norm), precision, r_norm);

    BFT_FREE(resr);
  }

  cs_timer_stats_switch(t_top_id);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_sles_t_tot, &t0, &t1);

  return state;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free sparse linear equation solver setup.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for
 * logging (especially performance data).
 *
 * \param[in, out]  sles  pointer to solver object
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_free(cs_sles_t  *sles)
{
  if (sles != NULL) {
    if (sles->free_func != NULL)
      sles->free_func(sles->context);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy the definition of a sparse linear equation solver to another.
 *
 * The intended use of this function is to allow associating different
 * solvers to related systems, so as to differentiate logging, while using
 * the same settings by default.
 *
 * If the source solver does not provide a \ref cs_sles_copy_t function,
 * No modification is done to the solver. If the copy function is available,
 * the context is copied, as are the matching function pointers.
 *
 * If previous settings have been defined and used, they are saved as
 * per \ref cs_sles_define.
 *
 * \param[in, out]  dest  pointer to destination solver object
 * \param[in]       src   pointer to source solver object
 *
 * \return  0 in case of success, 1 in case of failure
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_copy(cs_sles_t        *dest,
             const cs_sles_t  *src)
{
  int retval = 1;

  /* If no copy function is available or source does not have a
     context yet, we can do nothing */

  if (src->copy_func == NULL)
    return retval;

  /* Check if system was previously defined and used,
     and save info for future logging in this case */

  if (dest->context != NULL) {
    if (dest->n_calls > 0  && dest->log_func != NULL)
      _save_system_info(dest);
    else if (dest->destroy_func != NULL)
      dest->destroy_func(&(dest->context));
  }

  dest->type_id = src->type_id;
  dest->verbosity = src->verbosity;

  /* Now define options */
  dest->context = src->copy_func(src->context);
  dest->setup_func = src->setup_func;
  dest->solve_func = src->solve_func;
  dest->free_func = src->free_func;
  dest->log_func = src->log_func;
  dest->copy_func = src->copy_func;
  dest->destroy_func = src->destroy_func;

  if (dest->context != NULL)
    retval = 0;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a convergence error handler to a given sparse linear
 *        equation solver.
 *
 * The error will be called whenever convergence fails. To dissassociate
 * the error handler, this function may be called with \p handler = NULL.
 *
 * The association will only be successful if the matching solver
 * has already been defined.
 *
 * \param[in, out]  sles                pointer to solver object
 * \param[in]       error_handler_func  pointer to convergence error
 *                                      handler function
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_error_handler(cs_sles_t                *sles,
                          cs_sles_error_handler_t  *error_handler_func)
{
  if (sles != NULL)
    sles->error_func = error_handler_func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to default sparse linear solver definition function.
 *
 * The associated function will be used to provide a definition when
 * \ref cs_sles_setup or \ref cs_sles_solve is used for a system for which no
 * matching call to \ref cs_sles_define has been done.
 *
 * \return  define_func pointer to default definition function
 */
/*----------------------------------------------------------------------------*/

cs_sles_define_t  *
cs_sles_get_default_define(void)
{
  return _cs_sles_define_default;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set default sparse linear solver definition function.
 *
 * The provided function will be used to provide a definition when
 * \ref cs_sles_setup or \ref cs_sles_solve is used for a system for which no
 * matching call to \ref cs_sles_define has been done.
 *
 * \param[in]  define_func pointer to default definition function
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_default_define(cs_sles_define_t  *define_func)
{
  _cs_sles_define_default = define_func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set default verbosity definition function.
 *
 * The provided function will be used to define the verbosity when
 * \ref cs_sles_find_or_add is called.
 *
 * \param[in]  verbosity_func pointer to default verbosity function
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_default_verbosity(cs_sles_verbosity_t  *verbosity_func)
{
  _cs_sles_default_verbosity = verbosity_func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output default post-processing data for failed system convergence.
 *
 * \param[in]       name           variable name
 * \param[in]       mesh_id        id of error output mesh, or 0 if none
 * \param[in]       a              linear equation matrix
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             current system solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_def(const char          *name,
                              int                  mesh_id,
                              const cs_matrix_t   *a,
                              const cs_real_t     *rhs,
                              cs_real_t           *vx)
{
  if (mesh_id != 0) {

    const cs_mesh_t *mesh = cs_glob_mesh;

    char base_name[32], val_name[32];

    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a);
    const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);

    /* Check for mesh location */

    int location_id = 0, flag = 0;
    if (n_rows == mesh->n_cells)
      location_id = CS_MESH_LOCATION_CELLS;
    else if (n_rows == mesh->n_vertices)
      location_id = CS_MESH_LOCATION_VERTICES;
    flag = location_id;
    cs_parall_max(1, CS_INT_TYPE, &flag);
    if (flag == location_id)
      flag = 0;
    else
      flag =1;
    cs_parall_max(1, CS_INT_TYPE, &flag);
    if (flag != 0)
      return;

    /* Now generate output */

    cs_real_t *val;
    BFT_MALLOC(val, n_cols*diag_block_size, cs_real_t);

    for (int val_id = 0; val_id < 5; val_id++) {

      switch(val_id) {

      case 0:
        strcpy(base_name, "Diag");
        cs_matrix_copy_diagonal(a, val);
        break;

      case 1:
        strcpy(base_name, "RHS");
        memcpy(val, rhs, n_rows*diag_block_size*sizeof(cs_real_t));
        break;

      case 2:
        strcpy(base_name, "X");
        memcpy(val, vx, n_rows*diag_block_size*sizeof(cs_real_t));
        break;

      case 3:
        strcpy(base_name, "Residual");
        _residual(n_rows*diag_block_size,
                  a,
                  rhs,
                  vx,
                  val);
        break;

      case 4:
        strcpy(base_name, "Diag_Dom");
        cs_matrix_diag_dominance(a, val);
        break;

      }

      if (strlen(name) + strlen(base_name) < 31) {
        strcpy(val_name, base_name);
        strcat(val_name, "_");
        strcat(val_name, name);
      }
      else {
        strncpy(val_name, base_name, 31);
        val_name[31] = '\0';
      }

      cs_sles_post_output_var(val_name,
                              mesh_id,
                              location_id,
                              CS_POST_WRITER_ERRORS,
                              diag_block_size,
                              val);
    }

    BFT_FREE(val);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output post-processing variable related to system convergence.
 *
 * \param[in]       name             variable name
 * \param[in]       mesh_id          id of error output mesh, or 0 if none
 * \param[in]       location_id      mesh location id (cells or vertices)
 * \param[in]       writer_id        id of specified associated writer, or
 *                                   \ref CS_POST_WRITER_ALL_ASSOCIATED for all
 * \param[in]       diag_block_size  block size for diagonal
 * \param[in, out]  var              variable values
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_post_output_var(const char      *name,
                        int              mesh_id,
                        int              location_id,
                        int              writer_id,
                        int              diag_block_size,
                        cs_real_t        var[])
{
  if (mesh_id != 0) {

    int _diag_block_size[4] = {1, 1, 1, 1};

    const cs_mesh_t *mesh = cs_glob_mesh;
    const cs_time_step_t *ts = cs_glob_time_step;

    size_t n_non_norm;
    cs_lnum_t n_rows = 0;
    if (location_id == CS_MESH_LOCATION_CELLS)
      n_rows = mesh->n_cells;
    else if (location_id == CS_MESH_LOCATION_VERTICES)
      n_rows = mesh->n_vertices;

    cs_real_t *val_type;

    assert(_diag_block_size[0] == _diag_block_size[1]); /* no padding */

    if (diag_block_size > 1) {
      int i;
      for (i = 0; i < 3; i++) /* will become false if padding is used */
        _diag_block_size[i] = diag_block_size;
      _diag_block_size[3] = diag_block_size*diag_block_size;
    }

    BFT_MALLOC(val_type, _diag_block_size[1]*n_rows, cs_real_t);

    n_non_norm = _value_type(_diag_block_size[1]*n_rows, var, val_type);

    if (location_id == CS_MESH_LOCATION_CELLS)
      cs_post_write_var(mesh_id,
                        writer_id,
                        name,
                        _diag_block_size[0],
                        true, /* interlace */
                        true, /* use parents */
                        CS_POST_TYPE_cs_real_t,
                        var,
                        NULL,
                        NULL,
                        ts);
    else if (location_id == CS_MESH_LOCATION_VERTICES)
      cs_post_write_vertex_var(mesh_id,
                               writer_id,
                               name,
                               _diag_block_size[0],
                               true, /* interlace */
                               true, /* use parents */
                               CS_POST_TYPE_cs_real_t,
                               var,
                               ts);

    int flag = (n_non_norm > 0) ? 1 : 0;
    cs_parall_max(1, CS_INT_TYPE, &flag);

    if (flag > 0) {

      char type_name[32];
      size_t l = 31 - strlen("_fp_type");

      strncpy(type_name, name, l);
      type_name[31] = '\0';

      strcat(type_name, "_fp_type");

      if (location_id == CS_MESH_LOCATION_CELLS)
        cs_post_write_var(mesh_id,
                          writer_id,
                          type_name,
                          _diag_block_size[0],
                          true, /* interlace */
                          true, /* use parents */
                          CS_POST_TYPE_cs_real_t,
                          val_type,
                          NULL,
                          NULL,
                          ts);
      else if (location_id == CS_MESH_LOCATION_VERTICES)
        cs_post_write_vertex_var(mesh_id,
                                 writer_id,
                                 name,
                                 _diag_block_size[0],
                                 true, /* interlace */
                                 true, /* use parents */
                                 CS_POST_TYPE_cs_real_t,
                                 var,
                                 ts);

    }

    BFT_FREE(val_type);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return base name associated to a field id, name couple.
 *
 * This is simply a utility function which will return its name argument
 * if f_id < 0, and the associated field's name or label otherwise.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 *
 * \return  pointer to base name associated to the field id, name couple
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_base_name(int          f_id,
                  const char  *name)
{
  const char *sles_name = name;

  if (f_id > -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    sles_name = cs_field_get_label(f);
  }

  return sles_name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return name associated to a field id, name couple.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 *
 * \return  pointer to name associated to the field id, name couple
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_name(int          f_id,
             const char  *name)
{
  const cs_sles_t *sles = cs_sles_find_or_add(f_id, name);

  if (sles->name != NULL)
    return sles->name;
  else
    return cs_sles_base_name(f_id, name);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
