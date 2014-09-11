#ifndef __CS_SLES_H__
#define __CS_SLES_H__

/*============================================================================
 * Sparse Linear Equation Solver driver
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_log.h"
#include "cs_halo_perio.h"
#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convergence status
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SLES_DIVERGED = -3,
  CS_SLES_BREAKDOWN = -2,
  CS_SLES_MAX_ITERATION = -1,
  CS_SLES_ITERATING = 0,
  CS_SLES_CONVERGED = 1

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
 * residue/r_norm <= precision, residue being the L2 norm of a.vx-rhs.
 *
 * parameters:
 *   context       <-> pointer to solver context
 *   name          <-- pointer to name of linear system
 *   a             <-- matrix
 *   verbosity     <-- associated verbosity
 *   rotation_mode <-- halo update option for rotational periodicity
 *   precision     <-- solver precision
 *   r_norm        <-- residue normalization
 *   n_iter        --> number of "equivalent" iterations
 *   residue       --> residue
 *   rhs           <-- right hand side
 *   vx            <-- system solution
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
                   cs_halo_rotation_t   rotation_mode,
                   double               precision,
                   double               r_norm,
                   int                 *n_iter,
                   double              *residue,
                   const cs_real_t     *rhs,
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
 * An error handler may be associated with a given solver context using
 * cs_sles_set_error_handler(), in which case it will be called whenever
 * convergence fails.
 *
 * Remark: in the advent of re-trying with different parameters,
 * it is recommended that either the parameters be sufficiently similar that
 * performance logging will not be affected, so the info reported to
 * the user is not biased. Complex strategies involving different
 * solver types should not be based on the use of an error handler, but
 * built-into the solver function, with appropriate performance logging
 * (though they may use an error handler in case of final failure).
 *
 * parameters:
 *   context       <-> pointer to solver context
 *   state         <-- convergence state
 *   name          <-- name of linear system
 *   a             <-- matrix
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   rhs           <-- Right hand side
 *   vx            <-- System solution
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_error_handler_t) (void                         *context,
                           cs_sles_convergence_state_t   state,
                           const char                   *name,
                           const cs_matrix_t            *a,
                           cs_halo_rotation_t            rotation_mode,
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
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Log sparse linear equation solver info
 *
 * parameters:
 *   log_type <-- log type (CS_LOG_SETUP or CS_LOG_PERFORMANCE)
 *----------------------------------------------------------------------------*/

void
cs_sles_log(cs_log_t  log_type);

/*----------------------------------------------------------------------------
 * Return pointer to linear system object, based on matching field id or
 * system name.
 *
 * If this system did not previously exist, NULL is returned.
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   pointer to associated linear system object, or NULL
 *----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_find(int          f_id,
             const char  *name);

/*----------------------------------------------------------------------------
 * Return pointer to linear system object, based on matching field id or
 * system name.
 *
 * If this system did not previously exist, it is created and added to
 * to the list of "known" systems. In this case, it will be usable
 * only if cs_sles_define() is called for the same field id and name
 * (in which case calling the present function is redundant), or if
 * cs_sles_set_sefault_define() has been previously used to define
 * the default solver policy.
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   pointer to associated linear system object.
 *----------------------------------------------------------------------------*/

cs_sles_t *
cs_sles_find_or_add(int          f_id,
                    const char  *name);

/*----------------------------------------------------------------------------
 * Temporarily replace field id with name for matching calls
 * to cs_sles_setup(), cs_sles_solve(), cs_sles_free(), and other operations
 * involving access through a field id.
 *
 * Deprecated.
 *   This function is provided to allow some peculiar
 *   calling sequences, in which codits is called with a nonzero
 *   ivar value, but specific solver options must still be set.
 *   In the future, a more consistent mechanism (using a zero ivar
 *   value or designing a cleaner method to handle those exceptional cases)
 *   is preferred. As such, only a stack depth of 1 is allowed.
 *
 * parameters:
 *   f_id  associated field id, or < 0
 *   name  associated name if f_id < 0, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_push(int          f_id,
             const char  *name);

/*----------------------------------------------------------------------------
 * Restore behavior temporarily modified by cs_slesh_push().
 *
 * Deprecated.
 *   This function matches cs_sles_push(), which is deprecated.
 *
 * parameters:
 *   f_id  associated field id, or < 0
 *----------------------------------------------------------------------------*/

void
cs_sles_pop(int  f_id);

/*----------------------------------------------------------------------------
 * Define sparse linear equation solver for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. If it existed previously,
 *
 * The context pointer is used to point to a structure adapted to the function
 * pointers given here, and combined with those functions, allows using
 * both built-in, external, or user-defined solvers.
 *
 * It is recommended the context type name provided here directly relate
 * to the associated structure type (for example, "cs_sles_it_t" or
 * "cs_multigrid_t").
 *
 * parameters:
 *   f_id         <-- associated field id, or < 0
 *   name         <-- associated name if f_id < 0, or NULL
 *   context      <-> pointer to solver context management structure
 *                    (cs_sles subsystem becomes owner)
 *   setup_func   <-- pointer to system setup function
 *   solve_func   <-- pointer to system solution function
 *                    (also calls setup_func if not done yet
 *                    or free_func called since last solve)
 *   free_func    <-- pointer function freeing system setup
 *   log_func     <-- pointer to system info logging function
 *                    (optional, but recommended)
 *   copy_func    <-- pointer to settings copy function (optional)
 *   destroy_func <-- pointer to function destroying solver context
 *                                (called with cs_sles_finalize() or with
 *                                a new call to this function for the same
 *                                system)
 *
 * returns:
 *   pointer to associated linear system object
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Set the verbosity for a given linear equation solver.
 *
 * This verbosity will be used by cs_sles_setup() and cs_sles_solve().
 *
 * By default, the verbosity is set to 0, or the value returned by the
 * function set with cs_sles_set_default_define().
 *
 * parameters
 *   sles      <-> pointer to solver object
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sles_set_verbosity(cs_sles_t  *sles,
                      int         verbosity);

/*----------------------------------------------------------------------------
 * Return type name of solver context.
 *
 * The returned string is intended to help determine which type is associated
 * with the void * pointer returned by cs_sles_get_context() for a given
 * solver definition, so as to be able to call additional specific functions
 * beyond the generic functions assigned to a cs_sles_t object.
 *
 * If no type_name string was associated to the solver upon its definition by
 * cs_sles_define(), or it has not been defined yet, the string returned
 * is "<undefined>". It is recommended the type name provided
 * cs_sles_define() directly relate to the associated structure type
 * (for example, "cs_sles_it_t" or "cs_multigrid_t").
 *
 * parameters
 *   sles  pointer to solver object
 *
 * returns:
 *   pointer to linear system solver specific type name
 *----------------------------------------------------------------------------*/

const char *
cs_sles_get_type(cs_sles_t  *sles);

/*----------------------------------------------------------------------------
 * Return pointer to solver context structure pointer.
 *
 * The context structure depends on the type of solver used, which may in
 * turn be determined by the string returned by cs_sles_get_type().
 * If may be used by appropriate functions specific to that type.
 *
 * parameters
 *   sles  pointer to solver object
 *
 * returns:
 *   pointer to solver-specific linear system info and context
*----------------------------------------------------------------------------*/

void *
cs_sles_get_context(cs_sles_t  *sles);

/*----------------------------------------------------------------------------
 * Setup sparse linear equation solver.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * parameters:
 *   sles <-> pointer to solver object
 *   a    <-- matrix
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

void
cs_sles_setup(cs_sles_t          *sles,
              const cs_matrix_t  *a);

/*----------------------------------------------------------------------------
 * General sparse linear system resolution.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * Note that if cs_sles_setup() was previously called for this
 * system, and cs_sles_free() has not been called since, the matrix
 * provided should be the same. The optional separation between the
 * two stages is intended to allow amortizing the cost of setup
 * over multiple solutions.
 *
 * The system is considered to have converged when
 * residue/r_norm <= precision, residue being the L2 norm of a.vx-rhs.
 *
 * parameters:
 *   sles          <-> pointer to solver object
 *   a             <-- matrix
 *   rotation_mode <-- halo update option for rotational periodicity
 *   precision     <-- solver precision
 *   r_norm        <-- residue normalization
 *   n_iter        --> number of "equivalent" iterations
 *   residue       --> residue
 *   rhs           <-- right hand side
 *   vx            <-> system solution
 *   aux_size      <-- size of aux_vectors (in bytes)
 *   aux_vectors   --- optional working area (internal allocation if NULL)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve(cs_sles_t           *sles,
              const cs_matrix_t   *a,
              cs_halo_rotation_t   rotation_mode,
              double               precision,
              double               r_norm,
              int                 *n_iter,
              double              *residue,
              const cs_real_t     *rhs,
              cs_real_t           *vx,
              size_t               aux_size,
              void                *aux_vectors);

/*----------------------------------------------------------------------------
 * Free sparse linear equation solver setup.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for
 * logging (especially performance data).
 *
 * parameters:
 *   sles <-> pointer to solver object
 *----------------------------------------------------------------------------*/

void
cs_sles_free(cs_sles_t  *sles);

/*----------------------------------------------------------------------------
 * Copy the definition of a sparse linear equation solver to another.
 *
 * The intended use of this function is to allow associating different
 * solvers to related systems, so as to differentiate logging, while using
 * the same setttings by default.
 *
 * If the source solver does not provide a cs_sles_copy_t function,
 * No modification is done to the solver. If the copy function is available,
 * the context is copied, as are the matching function pointers.
 *
 * If previous settings have been defined and used, they are saved as
 * per cs_sles_define().
 *
 * parameters:
 *   dest <-> pointer to destination solver object
 *   src  <-- pointer to source solver object
 *
 * returns:
 *   0 in case of success, 1 in case of failure
 *----------------------------------------------------------------------------*/

int
cs_sles_copy(cs_sles_t        *dest,
             const cs_sles_t  *src);

/*----------------------------------------------------------------------------
 * Associate a convergence error handler to a given sparse linear
 * equation solver.
 *
 * The error will be called whenever convergence fails. To dissassociate
 * the error handler, this function may be called with \p handler = NULL.
 *
 * The association will only be successful if the matching solver
 * has already been defined.
 *
 * parameters:
 *   sles               <-> pointer to solver object
 *   error_handler_func <-- pointer to convergence error handler function
 *----------------------------------------------------------------------------*/

void
cs_sles_set_error_handler(cs_sles_t                *sles,
                          cs_sles_error_handler_t  *error_handler_func);

/*----------------------------------------------------------------------------
 * Set default sparse linear solver definition function.
 *
 * The provided function will be used to provide a definition when
 * cs_sles_setup() or cs_sles_solve() is used for a system for which no
 * matching call to cs_sles_define() has been done.
 *
 * parameters:
 *   define_func <-- pointer to default definition function
 *----------------------------------------------------------------------------*/

void
cs_sles_set_default_define(cs_sles_define_t  *define_func);

/*----------------------------------------------------------------------------
 * Set default verbosity definition function.
 *
 * The provided function will be used to define the verbosity when
 * cs_sles_default_verbosity() is called.
 *
 * parameters:
 *   verbosity_func <-- pointer to default verbosity function
 *----------------------------------------------------------------------------*/

void
cs_sles_set_default_verbosity(cs_sles_verbosity_t  *verbosity_func);

/*----------------------------------------------------------------------------
 * Test if a linear system needs solving or if the residue is already
 * within convergence criteria.
 *
 * parameters:
 *   solver_name <-- name of the solver calling the test
 *   system_name <-- name of the linear system tested
 *   verbosity   <-- verbosity level
 *   precision   <-- solver precision
 *   r_norm      <-- residue normalization
 *   residue     <-- residue
 *
 * returns:
 *   1 if solving is required, 0 if the residue is already zero within
 *   tolerance criteria (precision of residue normalization)
 *----------------------------------------------------------------------------*/

int
cs_sles_needs_solving(const char  *solver_name,
                      const char  *system_name,
                      int          verbosity,
                      double       precision,
                      double       r_norm,
                      double       residue);

/*----------------------------------------------------------------------------
 * Output default post-processing data for failed system convergence.
 *
 * parameters:
 *   name             <-- variable name
 *   mesh_id          <-- id of error output mesh, or 0 if none
 *   rotation_mode    <-- halo update option for rotational periodicity
 *   a                <-- linear equation matrix
 *   rhs              <-- right hand side
 *   vx               <-> current system solution
 *----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_def(const char          *name,
                              int                  mesh_id,
                              cs_halo_rotation_t   rotation_mode,
                              const cs_matrix_t   *a,
                              const cs_real_t     *rhs,
                              cs_real_t           *vx);

/*----------------------------------------------------------------------------
 * Output post-processing variable for failed system convergence.
 *
 * parameters:
 *   name            <-- variable name
 *   diag_block_size <-- block size for diagonal
 *   mesh_id         <-- id of error output mesh, or 0 if none
 *   var             <-> variable values
 *----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_var(const char  *name,
                              int          mesh_id,
                              int          diag_block_size,
                              cs_real_t   *var);

/*----------------------------------------------------------------------------*
 * Return base name associated to a field id, name couple.
 *
 * This is simply a utility function which will return its name argument
 * if f_id < 0, and the associated field's name or label otherwise.
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   pointer to associated linear system object, or NULL
 *----------------------------------------------------------------------------*/

const char *
cs_sles_base_name(int          f_id,
                  const char  *name);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_H__ */
