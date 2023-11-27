#ifndef __CS_SLES_PC_H__
#define __CS_SLES_PC_H__

/*============================================================================
 * Sparse Linear Equation Solver Preconditioner driver
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

  CS_SLES_PC_DIVERGED = -2,
  CS_SLES_PC_BREAKDOWN = -1,
  CS_SLES_PC_MAX_ITERATION = 0,
  CS_SLES_PC_CONVERGED = 1

} cs_sles_pc_state_t;

/* General linear solver context (opaque) */

typedef struct _cs_sles_pc_t  cs_sles_pc_t;

/*----------------------------------------------------------------------------
 * Function pointer returning the type name of a preconditioner context.
 *
 * The context structure depends on the type of preconditioner used,
 * which may in turn be determined by the string returned by
 * cs_sles_pc_get_type() and cs_sles_pc_get_type_name().
 * If may be used by appropriate functions specific to that type.
 *
 * parameters:
 *   context   <-- pointer to preconditioner-specific context
 *   logging   <-- if true, a name appropritate to logging
 *                 (possibly translated) is returned; if false,
 *                 a canonical name is returned.
 *----------------------------------------------------------------------------*/

typedef const char *
(cs_sles_pc_get_type_t) (const void  *context,
                         bool         logging);

/*----------------------------------------------------------------------------
 * Function pointer for pre-resolution setup of a preconditioner context.
 *
 * This setup may include building a multigrid hierarchy, for example.
 *
 * parameters:
 *   context   <-> pointer to preconditioner context
 *   name      <-- pointer to name of associated linear system
 *   a         <-- matrix
 *   accel     <-- use accelerator version ?
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_pc_setup_t) (void               *context,
                      const char         *name,
                      const cs_matrix_t  *a,
                      bool                accel,
                      int                 verbosity);

/*----------------------------------------------------------------------------
 * Function pointer for setting of the required tolerance for preconditioners
 * involving an iterative solver.
 *
 * This will usually not be relevant to non-iterative preconditioners,
 * for which this type of function does not need to be defined.
 *
 * The preconditioner is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   precision     <-- preconditioner precision
 *   r_norm        <-- residual normalization
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_pc_tolerance_t) (void    *context,
                          double   precision,
                          double   r_norm);

/*----------------------------------------------------------------------------
 * Function pointer for application of a preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

typedef cs_sles_pc_state_t
(cs_sles_pc_apply_t) (void                *context,
                      const cs_real_t     *x_in,
                      cs_real_t           *x_out);

/*----------------------------------------------------------------------------
 * Function pointer for freeing of a preconditioner's context data.
 *
 * Note that this function should free resolution-related data, such as
 * multigrid hierarchy and any other temporary arrays or
 * objects required for application, but should not free the whole context,
 * as info used for logging (especially performance data) should be
 * maintained.
 *
 * parameters:
 *   context <-> pointer to preconditioner context
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_pc_free_t) (void  *context);

/*----------------------------------------------------------------------------
 * Function pointer for logging of linear preconditioner setup,
 * history and performance data.
 *
 * This function will indirectly  be called for each preconditioner when
 * cs_sles_finalize() is called.
 *
 * parameters:
 *   context  <-- pointer to preconditioner context
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_pc_log_t) (const void  *context,
                    cs_log_t     log_type);

/*----------------------------------------------------------------------------
 * Function pointer for creation of a preconditioner context based on the
 * copy of another.
 *
 * The new context copies the settings of the copied context, but not
 * its setup data and logged info, such as performance data.
 *
 * This type of function is optional, but enables associating different
 * preconditioners to related systems (to differentiate logging) while using
 * the same settings by default.
 *
 * parameters:
 *   context  <-- context to clone
 *
 * returns:
 *   pointer to newly created context
 *----------------------------------------------------------------------------*/

typedef void *
(cs_sles_pc_clone_t) (const void  *context);

/*----------------------------------------------------------------------------
 * Function pointer for destruction of a preconditioner context.
 *
 * This function should free all context data.
 *
 * parameters:
 *   context <-> pointer to preconditioner context
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_pc_destroy_t) (void  **context);

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
 * \brief Log preconditioner setup, history and performance data.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for
 * logging (especially performance data).
 *
 * \param[in, out]  pc        pointer to preconditioner object
 * \param[in]       log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_log(cs_sles_pc_t  *pc,
               cs_log_t       log_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define sparse linear equation preconditioner.
 *
 * The context pointer is used to point to a structure adapted to the function
 * pointers given here, and combined with those functions, allows using
 * both built-in, external, or user-defined preconditioners.
 *
 * \param[in, out]  context         pointer to preconditioner context structure
 *                                  (cs_sles_pc subsystem becomes owner)
 * \param[in]       get_type_func   pointer to function returning type name
 * \param[in]       setup_func      pointer to preconditioner setup function
 * \param[in]       tolerance_func  pointer to tolerance setting functio
 * \param[in]       apply_func      pointer to preconditioner application
 *                                  function (also calls setup_func if not done
 *                                  yet or free_func called since last apply)
 * \param[in]       free_func       pointer function freeing system setup
 * \param[in]       log_func        pointer to system info logging function
                                    (optional, but recommended)
 * \param[in]       clone_func      pointer to settings clone function
 * \param[in]       destroy_func    pointer to function destroying
 *                                  preconditioner context
 *
 * \return  pointer to associated preconditioner object
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_define(void                    *context,
                  cs_sles_pc_get_type_t   *get_type_func,
                  cs_sles_pc_setup_t      *setup_func,
                  cs_sles_pc_tolerance_t  *tolerance_func,
                  cs_sles_pc_apply_t      *apply_func,
                  cs_sles_pc_free_t       *free_func,
                  cs_sles_pc_log_t        *log_func,
                  cs_sles_pc_clone_t      *clone_func,
                  cs_sles_pc_destroy_t    *destroy_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a sparse linear equation preconditioner.
 *
 *
 * \param[in, out]  pc   pointer to preconditioner context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_destroy(cs_sles_pc_t **pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new preconditioner context based on the copy of another.
 *
 * The intended use of this function is to allow associating different
 * preconditioners to related systems, so as to allow simultaneous setups
 * and differentiate logging, while using the same settings by default.
 *
 * If no preconditioner (i.e. NULL) is passed, it will return NULL.
 *
 * \param[in]       src   pointer to source preconditioner object
 *
 * \return  pointer to new preconditioner object, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_clone(const cs_sles_pc_t  *src);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return type name of preconditioner context.
 *
 * The returned string is intended to help determine which type is associated
 * with the void * pointer returned by \ref cs_sles_pc_get_context for a given
 * preconditioner definition, so as to be able to call additional specific
 * functions beyond the generic functions assigned to a cs_sles_pc_t object.
 *
 * \param[in]  pc  pointer to preconditioner object
 *
 * \return  pointer to linear system preconditioner specific type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_pc_get_type(cs_sles_pc_t  *pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return type name of preconditioner context.
 *
 * The returned string is intended mainly for logging.
 *
 * \param[in]  pc  pointer to preconditioner object
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_pc_get_type_name(cs_sles_pc_t  *pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to preconditioner context structure pointer.
 *
 * The context structure depends on the type of preconditioner used, which may
 * in turn be determined by the string returned by cs_sles_pc_get_type().
 * If may be used by appropriate functions specific to that type.
 *
 * \param[in]  pc  pointer to preconditioner object
 *
 * \return  pointer to preconditioner-specific info and context
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_pc_get_context(cs_sles_pc_t  *pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the function used to apply a preconditioner.
 *
 * This allows calling the preconditioner with one less level of indirection.
 *
 * \param[in]  pc  pointer to preconditioner object
 *
 * \return  preconditioner apply function
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_apply_t *
cs_sles_pc_get_apply_func(const cs_sles_pc_t *pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the required tolerance for preconditioners involving an
 *        iterative solver.
 *
 * This will usually not be relevant to non-iterative preconditioners,
 * in which case this is a no-op.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * The system is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 *
 * \param[in, out]  pc             pointer to preconditioner object
 * \param[in]       precision      preconditioner precision
 * \param[in]       r_norm         residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_set_tolerance(cs_sles_pc_t  *pc,
                         double         precision,
                         double         r_norm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup sparse linear equation preconditioner.
 *
 * Use of this function is optional: if a \ref cs_sles_solve is called
 * for the same system before this function is called, the latter will be
 * called automatically.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * \param[in, out]  pc         pointer to preconditioner object
 * \param[in]       name       linear system name
 * \param[in]       a          matrix
 * \param[in]       accel      use accelerator version ?
 * \param[in]       verbosity  verbosity level
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_setup(cs_sles_pc_t       *pc,
                 const char         *name,
                 const cs_matrix_t  *a,
                 bool                accel,
                 int                 verbosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a preconditioner.
 *
 * If no options were previously provided for the matching system,
 * default options will be used.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * \param[in, out]  pc             pointer to preconditioner object
 * \param[in]       x_in           input vector
 * \param[in, out]  x_out          input/output vector
 *
 * \return  preconditioner application status
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_apply(cs_sles_pc_t        *pc,
                 cs_real_t           *x_in,
                 cs_real_t           *x_out);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free preconditioner setup.
 *
 * This function frees resolution-related data, such as multigrid hierarchy,
 * preconditioning, and any other temporary arrays or objects required for
 * resolution, but maintains context information such as that used for
 * logging (especially performance data).
 *
 * \param[in, out]  pc  pointer to preconditioner object
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_pc_free(cs_sles_pc_t  *pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an "identity" (or null) preconditioner.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_none_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Jacobi preconditioner.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_jacobi_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a polynomial preconditioner of degree 1.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_poly_1_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a polynomial preconditioner of degree 2.
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_sles_pc_poly_2_create(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_PC_H__ */
