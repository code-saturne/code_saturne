#ifndef CS_SLES_H
#define CS_SLES_H

/*============================================================================
 * Sparse Linear Equation Solver driver
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_log.h"
#include "base/cs_halo_perio.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convergence status
 *----------------------------------------------------------------------------*/

/*! Convergence status indicator. */
typedef enum {

  CS_SLES_DIVERGED = -3,      /*!< The solver has diverged */
  CS_SLES_BREAKDOWN = -2,     /*!< The solver has broken down, and cannot make
                                   any more progress*/
  CS_SLES_MAX_ITERATION = -1, /*!< Maximum number of iterations has been
                                   reached, without reaching the required
                                   convergence */
  CS_SLES_ITERATING = 0,      /*!< The solver is iterating */
  CS_SLES_CONVERGED = 1       /*!< The solver has converged */

} cs_sles_convergence_state_t;

/* General linear solver context (opaque) */

typedef struct _cs_sles_t  cs_sles_t;

/*----------------------------------------------------------------------------
 * Function pointer for pre-resolution setup of a linear system solvers's
 * context.
 *
 * This setup may include building a multigrid hierarchy, or a preconditioner.
 *
 * Use of this type of function is optional: the context is expected to
 * maintain state, so that if a cs_sles_solve_t function is called before a
 * cs_sles_setup_t function, the latter will be called automatically.
 *
 * parameters:
 *   context   <-> pointer to solver context
 *   name      <-- pointer to name of linear system
 *   a         <-- matrix
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_setup_t) (void               *context,
                   const char         *name,
                   const cs_matrix_t  *a,
                   int                 verbosity);

/*----------------------------------------------------------------------------
 * Function pointer for resolution of a linear system.
 *
 * If the associated cs_sles_setup_t function has not been called before
 * this function, it will be called automatically.
 *
 * The solution context setup by this call (or that of the matching setup
 * function) will be maintained until the matching cs_sles_free_t function
 * is called.
 *
 * The matrix is not expected to change between successive calls, although
 * the right hand side may. If the matrix changes, the associated
 * cs_sles_setup_t or cs_sles_free_t function must be called between
 * solves.
 *
 * The system is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 *
 * parameters:
 *   context       <-> pointer to solver context
 *   name          <-- pointer to name of linear system
 *   a             <-- matrix
 *   verbosity     <-- associated verbosity
 *   precision     <-- solver precision
 *   r_norm        <-- residual normalization
 *   n_iter        --> number of "equivalent" iterations
 *   residual      --> residual
 *   rhs           <-- right hand side
 *   vx_ini        <-- initial solution (vx if nonzero, nullptr if zero)
 *   vx            --> system solution
 *   aux_size      <-- number of elements in aux_vectors
 *   aux_vectors   <-- optional working area (internal allocation if NULL)
 *
 * returns:
 *   convergence status
 *----------------------------------------------------------------------------*/

typedef cs_sles_convergence_state_t
(cs_sles_solve_t) (void                *context,
                   const char          *name,
                   const cs_matrix_t   *a,
                   int                  verbosity,
                   double               precision,
                   double               r_norm,
                   int                 *n_iter,
                   double              *residual,
                   const cs_real_t     *rhs,
                   cs_real_t           *vx_ini,
                   cs_real_t           *vx,
                   size_t               aux_size,
                   void                *aux_vectors);

/*----------------------------------------------------------------------------
 * Function pointer for freeing of a linear system's context data.
 *
 * Note that this function should free resolution-related data, such as
 * multigrid hierarchy, preconditioning, and any other temporary arrays or
 * objects required for resolution, but should not free the whole context,
 * as info used for logging (especially performance data) should be
 * maintained.
 *
 * parameters:
 *   context <-> pointer to solver context
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_free_t) (void  *context);

/*----------------------------------------------------------------------------
 * Function pointer for logging of linear solver setup,
 * history and performance data.
 *
 * This function will be called for each solver when cs_sles_finalize()
 * is called.
 *
 * parameters:
 *   context  <-- pointer to solver context
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_log_t) (const void  *context,
                 cs_log_t     log_type);

/*----------------------------------------------------------------------------
 * Function pointer for creation of a solver context based on the copy
 * of another.
 *
 * The new context copies the settings of the copied context, but not
 * its setup data and logged info, such as performance data.
 *
 * This type of function is optional, but enables associating different
 * solvers to related systems (to differentiate logging) while using
 * the same settings by default.
 *
 * parameters:
 *   context  <-- source context
 *
 * returns:
 *   pointer to newly created context
 *----------------------------------------------------------------------------*/

typedef void *
(cs_sles_copy_t) (const void  *context);

/*----------------------------------------------------------------------------
 * Function pointer for destruction of a linear system solver context.
 *
 * This function should free all context data, and will be called for each
 * system when cs_sles_finalize() is called.
 *
 * parameters:
 *   context <-> pointer to solver context
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_destroy_t) (void  **context);

/*----------------------------------------------------------------------------
 * Function pointer for handling of non-convegence when solving
 * a linear system.
 *
 * Such a function is optional, and may be used for a variety of purposes,
 * such as logging, postprocessing, re-trying with different parameters,
 * aborting the run, or any combination thereof.
 *
 * An error handler may be associated with a given solver using
 * cs_sles_set_error_handler(), in which case it will be called whenever
 * convergence fails.
 *
 * parameters:
 *   sles          <-> pointer to solver object
 *   state         <-- convergence status
 *   a             <-- matrix
 *   rhs           <-- Right hand side
 *   vx            <-- System solution
 *
 * returns:
 *   true if solve should be re-executed, false otherwise
 *----------------------------------------------------------------------------*/

typedef bool
(cs_sles_error_handler_t) (cs_sles_t                    *sles,
                           cs_sles_convergence_state_t   state,
                           const cs_matrix_t            *a,
                           const cs_real_t              *rhs,
                           cs_real_t                    *vx);

/*----------------------------------------------------------------------------
 * Function pointer for the default definition of a sparse
 * linear equation solver
 *
 * The function may be associated using cs_sles_set_default_define(), so
 * that it may provide a definition that will be used when
 * cs_sles_setup() or cs_sles_solve() is used for a system for which
 * no matching call to cs_sles_define() has been done.
 *
 * The function should call cs_sles_define() with arguments f_id
 * and name, and appropriately chosen function pointers.
 *
 * A pointer to the matrix of the system to be solved is also provided,
 * so that the corresponding information may be used to better choose
 * defaults.
 *
 * parameters:
 *   f_id          <-- associated field id, or < 0
 *   name          <-- associated name if f_id < 0, or NULL
 *   a             <-- Matrix
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_define_t) (int                 f_id,
                    const char         *name,
                    const cs_matrix_t  *a);

/*----------------------------------------------------------------------------
 * Function pointer for the default definition of a sparse
 * linear equation solver's verbosity
 *
 * The function may be associated using cs_sles_set_default_verbosity(), so
 * that it may provide a definition that will be used when
 * cs_sles_default_verbosity() is called.
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   default verbosity value
 *----------------------------------------------------------------------------*/

typedef int
(cs_sles_verbosity_t) (int          f_id,
                       const char  *name);

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Set the threshold value used in the detection of immediate exit
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_epzero(double  new_value);

/*----------------------------------------------------------------------------*/
/*
 * Get the current threshold value used in the detection of immediate exit
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_get_epzero(void);

/*----------------------------------------------------------------------------
 * Initialize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_finalize(void);

/*----------------------------------------------------------------------------*/
/*
 * Log sparse linear equation solver info
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_log(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/
/*
 * Return pointer to linear system object, based on matching field id or system
 * name.
 *
 * If this system did not previously exist, nullptr is returned.
 */
/*----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_find(int          f_id,
             const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Return pointer to linear system object, based on matching field id or system
 * name.
 *
 * If this system did not previously exist, it is created and added to to the
 * list of "known" systems. In this case, it will be usable only if
 * cs_sles_define() is called for the same field id and name (in which case
 * calling the present function is redundant), or if
 * cs_sles_set_sefault_define() has been previously used to define the default
 * solver policy.
 */
/*----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_find_or_add(int          f_id,
                    const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Temporarily replace field id with name for matching calls to cs_sles_setup,
 * cs_sles_solve, cs_sles_free, and other operations involving access through a
 * field id.
 *
 * This function is provided to allow some peculiar calling sequences, in which
 * cs_equation_iterative_solve_scalar is called with a given field id, but
 * specific solver options must still be set.  In the future, a cleaner method
 * to handle those exceptional cases would be preferred. As such, only a stack
 * depth of 1 is allowed.
 */
/*----------------------------------------------------------------------------*/

[[deprecated]] void
cs_sles_push(int          f_id,
             const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Restore behavior temporarily modified by cs_sles_push.
 *
 * This function matches cs_sles_push, which is deprecated.
 */
/*----------------------------------------------------------------------------*/

[[deprecated]] void
cs_sles_pop(int  f_id);

/*----------------------------------------------------------------------------*/
/*
 * Define sparse linear equation solver for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of "known"
 * systems.
 *
 * The context pointer is used to point to a structure adapted to the function
 * pointers given here, and combined with those functions, allows using both
 * built-in, external, or user-defined solvers.
 *
 * It is recommended the context type name provided here directly relate to the
 * associated structure type (for example, "cs_sles_it_t" or "cs_multigrid_t").
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
               cs_sles_destroy_t  *destroy_func);

/*----------------------------------------------------------------------------*/
/*
 * Set the verbosity for a given linear equation solver.
 *
 * This verbosity will be used by cs_sles_setup and cs_sles_solve.
 *
 * By default, the verbosity is set to 0, or the value returned by the function
 * set with cs_sles_set_default_define().
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_verbosity(cs_sles_t  *sles,
                      int         verbosity);

/*----------------------------------------------------------------------------*/
/*
 * Get the verbosity for a given linear equation solver.
 *
 * This verbosity will be used by cs_sles_setup and cs_sles_solve.
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_verbosity(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Activate postprocessing output for a given linear equation solver.
 *
 * This allows the output of the residual at the end of each solution series,
 * using a single postprocessing writer. By default, no output is activated.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_post_output(cs_sles_t  *sles,
                        int         writer_id);

/*----------------------------------------------------------------------------*/
/*
 * Return the id of the associated writer if postprocessing output is active
 * for a given linear equation solver.
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_post_output(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Return type name of solver context.
 *
 * The returned string is intended to help determine which type is associated
 * with the void * pointer returned by cs_sles_get_context for a given
 * solver definition, so as to be able to call additional specific functions
 * beyond the generic functions assigned to a cs_sles_t object.
 *
 * If no type_name string was associated to the solver upon its definition by
 * cs_sles_define, or it has not been defined yet, the string returned
 * is "<undefined>". It is recommended the type name provided
 * cs_sles_define directly relate to the associated structure type
 * (for example, "cs_sles_it_t" or "cs_multigrid_t").
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_get_type(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Return pointer to solver context structure pointer.
 *
 * The context structure depends on the type of solver used, which may in turn
 * be determined by the string returned by cs_sles_get_type().  If may be used
 * by appropriate functions specific to that type.
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_get_context(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Return field id associated with a given sparse linear equation solver.
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_f_id(const cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Return name associated with a given sparse linear equation solver.
 *
 * This is simply a utility function which will return its name argument if
 * f_id < 0, and the associated field's name otherwise.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_get_name(const cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Query if immediate_return ("no-op") is allowed when initial guess is zero
 * (solve by increments) and the RHS is already zero within the normalized
 * tolerance criteria.
 */
/*----------------------------------------------------------------------------*/

bool
cs_sles_get_allow_no_op(const cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Indicate if immediate_return ("no-op") is allowed when initial guess is zero
 * (solve by increments) and the RHS is already zero within the normalized
 * tolerance criteria.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_allow_no_op(cs_sles_t  *sles,
                        bool        allow_no_op);

/*----------------------------------------------------------------------------*/
/*
 * Setup sparse linear equation solver.
 *
 * Use of this function is optional: if a cs_sles_solve is called for the same
 * system before this function is called, the latter will be called
 * automatically.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup(cs_sles_t          *sles,
              const cs_matrix_t  *a);

/*----------------------------------------------------------------------------*/
/*
 * General sparse linear system resolution.
 *
 * If no options were previously provided for the matching system, default
 * options will be used.
 *
 * Note that if cs_sles_setup was previously called for this system, and
 * cs_sles_free has not been called since, the matrix provided should be the
 * same. The optional separation between the two stages is intended to allow
 * amortizing the cost of setup over multiple solutions.
 *
 * If the initial solution iz zero, setting vx_ini_0 to true may allow extra
 * optimizations with some solvers, and avoids requiring an upstream
 * initialization.
 *
 * The system is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve(cs_sles_t           *sles,
              const cs_matrix_t   *a,
              double               precision,
              double               r_norm,
              bool                 vx_ini_0,
              int                 *n_iter,
              double              *residual,
              const cs_real_t     *rhs,
              cs_real_t           *vx,
              size_t               aux_size,
              void                *aux_vectors);

/*----------------------------------------------------------------------------*/
/*
 * Free sparse linear equation solver setup.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for logging
 * (especially performance data).
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_free(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * Copy the definition of a sparse linear equation solver to another.
 *
 * The intended use of this function is to allow associating different solvers
 * to related systems, so as to differentiate logging, while using the same
 * settings by default.
 *
 * If the source solver does not provide a cs_sles_copy_t function, No
 * modification is done to the solver. If the copy function is available, the
 * context is copied, as are the matching function pointers.
 *
 * If previous settings have been defined and used, they are saved as per
 * cs_sles_define.
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_copy(cs_sles_t        *dest,
             const cs_sles_t  *src);

/*----------------------------------------------------------------------------*/
/*
 * Associate a convergence error handler to a given sparse linear equation
 * solver.
 *
 * The error will be called whenever convergence fails. To dissassociate the
 * error handler, this function may be called with \p handler = NULL.
 *
 * The association will only be successful if the matching solver
 * has already been defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_error_handler(cs_sles_t                *sles,
                          cs_sles_error_handler_t  *error_handler_func);

/*----------------------------------------------------------------------------*/
/*
 * Return pointer to default sparse linear solver definition function.
 *
 * The associated function will be used to provide a definition when
 * cs_sles_setup or cs_sles_solve is used for a system for which no matching
 * call to cs_sles_define has been done.
 */
/*----------------------------------------------------------------------------*/

cs_sles_define_t  *
cs_sles_get_default_define(void);

/*----------------------------------------------------------------------------*/
/*
 * Set default sparse linear solver definition function.
 *
 * The provided function will be used to provide a definition when
 * cs_sles_setup or cs_sles_solve is used for a system for which no matching
 * call to cs_sles_define has been done.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_default_define(cs_sles_define_t  *define_func);

/*----------------------------------------------------------------------------*/
/*
 * Set default verbosity definition function.
 *
 * The provided function will be used to define the verbosity when
 * cs_sles_find_or_add is called.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_default_verbosity(cs_sles_verbosity_t  *verbosity_func);

/*----------------------------------------------------------------------------*/
/*
 * Output default post-processing data for failed system convergence.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_def(const char          *name,
                              int                  mesh_id,
                              const cs_matrix_t   *a,
                              const cs_real_t     *rhs,
                              cs_real_t           *vx);

/*----------------------------------------------------------------------------*/
/*
 * Output post-processing variable related to system convergence.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_post_output_var(const char      *name,
                        int              mesh_id,
                        int              location_id,
                        int              writer_id,
                        int              diag_block_size,
                        cs_real_t        var[]);

/*----------------------------------------------------------------------------*/
/*
 * Return base name associated to a field id, name couple.
 *
 * This is simply a utility function which will return its name argument if
 * f_id < 0, and the associated field's name otherwise.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_base_name(int          f_id,
                  const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Return name associated to a field id, name couple.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_name(int          f_id,
             const char  *name);

/*----------------------------------------------------------------------------*/

#endif /* CS_SLES_H */
