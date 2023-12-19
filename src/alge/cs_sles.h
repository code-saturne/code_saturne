#ifndef __CS_SLES_H__
#define __CS_SLES_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_log.h"
#include "cs_halo_perio.h"
#include "cs_matrix.h"
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
                   double               precision,
                   double               r_norm,
                   int                 *n_iter,
                   double              *residual,
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
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the threshold value used in the detection of immediate exit
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_epzero(double  new_value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the current threshold value used in the detection of immediate
 *        exit
 *
 * \return the value of the threshold
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_get_epzero(void);

/*----------------------------------------------------------------------------
 * \brief Initialize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_initialize(void);

/*----------------------------------------------------------------------------
 * \brief Finalize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info
 *
 * \param[in]  log_type  log type (CS_LOG_SETUP or CS_LOG_PERFORMANCE)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_log(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/
/*
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
             const char  *name);

/*----------------------------------------------------------------------------*/
/*
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
                    const char  *name);

/*----------------------------------------------------------------------------*/
/*
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
             const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Restore behavior temporarily modified by \ref cs_sles_push.
 *
 * \deprecated This function matches \ref cs_sles_push, which is deprecated.
 *
 * \param[in]  f_id  associated field id, or < 0
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pop(int  f_id);

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
               cs_sles_destroy_t  *destroy_func);

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
                      int         verbosity);

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
cs_sles_get_verbosity(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
                        int         writer_id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return the id of the associated writer if postprocessing output
 *        is active for a given linear equation solver.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  id od associated writer, or 0
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_post_output(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_get_type(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_get_context(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return field id associated with a given sparse linear equation solver.
 *
 * \param[in]  sles  pointer to solver object
 *
 * \return  associated field id (or -1 if defined by name)
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_get_f_id(const cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_get_name(const cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_get_allow_no_op(const cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
                        bool        allow_no_op);

/*----------------------------------------------------------------------------*/
/*
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
              const cs_matrix_t  *a);

/*----------------------------------------------------------------------------*/
/*
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
              void                *aux_vectors);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_free(cs_sles_t  *sles);

/*----------------------------------------------------------------------------*/
/*
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
             const cs_sles_t  *src);

/*----------------------------------------------------------------------------*/
/*
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
                          cs_sles_error_handler_t  *error_handler_func);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_get_default_define(void);

/*----------------------------------------------------------------------------*/
/*
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
cs_sles_set_default_define(cs_sles_define_t  *define_func);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set default verbosity definition function.
 *
 * The provided function will be used to define the verbosity when
 * \ref cs_sles_find_or_add is called.
 *
 * \param[in]  verbosity_func pointer to default verbosity function
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_set_default_verbosity(cs_sles_verbosity_t  *verbosity_func);

/*----------------------------------------------------------------------------*/
/*
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
                              cs_real_t           *vx);

/*----------------------------------------------------------------------------*/
/*
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
                        cs_real_t        var[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return base name associated to a field id, name couple.
 *
 * This is simply a utility function which will return its name argument
 * if f_id < 0, and the associated field's name or label otherwise.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 *
 * \return  pointer to associated linear system object, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_base_name(int          f_id,
                  const char  *name);

/*----------------------------------------------------------------------------*/
/*
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
             const char  *name);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_H__ */
