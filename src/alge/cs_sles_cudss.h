#ifndef CS_SLES_CUDSS_H
#define CS_SLES_CUDSS_H

/*============================================================================
 * Sparse Linear Equation Solvers using cuDSS
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_base.h"
#include "alge/cs_matrix.h"
#include "alge/cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * cuDSS wrapper option flags
 */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* cuDSS linear solver context (opaque) */

typedef struct _cs_sles_cudss_t  cs_sles_cudss_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define and associate an cuDSS linear system solver
 * for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * cs_sles_define() and cs_sles_cudss_create().
 *
 * Note that this function returns a pointer directly to the cuDSS solver
 * management structure. This may be used to set further options.
 * If needed, cs_sles_find() may be used to obtain a pointer to the matching
 * cs_sles_t container.
 *
 * parameters:
 *   f_id         <-- associated field id, or < 0
 *   name         <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   pointer to newly created cuDSS solver info object.
 *----------------------------------------------------------------------------*/

cs_sles_cudss_t *
cs_sles_cudss_define(int          f_id,
                     const char  *name);

/*----------------------------------------------------------------------------
 * Create cuDSS linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * returns:
 *   pointer to newly created solver info object.
 *----------------------------------------------------------------------------*/

cs_sles_cudss_t *
cs_sles_cudss_create(void);

/*----------------------------------------------------------------------------
 * Create cuDSS linear system solver info and context
 * based on existing info and context.
 *
 * Most configuration parameters will be copied from the existing context,
 * though not all, since cuDSS does not provide a comprehensive way to
 * do this.
 *
 * parameters:
 *   context <-- pointer to reference info and context
 *               (actual type: cs_sles_cudss_t  *)
 *
 * returns:
 *   pointer to newly created solver info object
 *   (actual type: cs_sles_cudss_t  *)
 *----------------------------------------------------------------------------*/

void *
cs_sles_cudss_copy(const void  *context);

/*----------------------------------------------------------------------------
 * Destroy cuDSS linear system solver info and context.
 *
 * parameters:
 *   context  <-> pointer to cuDSS linear solver info
 *                (actual type: cs_sles_cudss_t  **)
 *----------------------------------------------------------------------------*/

void
cs_sles_cudss_destroy(void  **context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define additional cuDSS solver usage flags
 *
 * By default, the device will be used, but by calling this function
 * with "use_device = false", only the host will be used.
 *
 * \param[in, out]  context   pointer to cuDSS solver info and context
 * \param[in]       flags     flags (sum/bitwise of) for cuDSS usage options.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_set_flags(void  *context,
                        int    flags);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query additional cuDSS solver usage flags.
 *
 * \param[in]  context  pointer to cuDSS solver info and context
 *
 * \return  associated flags
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_cudss_get_flags(void  *context);

/*----------------------------------------------------------------------------
 * Setup cuDSS linear equation solver.
 *
 * parameters:
 *   context   <-> pointer to cuDSS linear solver info
 *                 (actual type: cs_sles_cudss_t  *)
 *   name      <-- pointer to system name
 *   a         <-- associated matrix
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sles_cudss_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity);

/*----------------------------------------------------------------------------
 * Call cuDSS linear equation solver.
 *
 * \warn The precision, r_norm, and n_iter parameters are ignored here.
 *       the matching configuration options should be set earlier, using
 *       the \ref cs_sles_cudss_set_config function
 *
 *
 * parameters:
 *   context       <-> pointer to cuDSS linear solver info
 *                     (actual type: cs_sles_cudss_t  *)
 *   name          <-- pointer to system name
 *   a             <-- matrix
 *   verbosity     <-- verbosity level
 *   precision     <-- solver precision
 *   r_norm        <-- residual normalization
 *   n_iter        --> number of iterations
 *   residual      --> residual
 *   rhs           <-- right hand side
 *   vx_ini        <-- initial system solution
 *                     (vx if nonzero, nullptr if zero)
 *   vx            <-> system solution
 *   aux_size      <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors   --- optional working area (internal allocation if NULL)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_cudss_solve(void                *context,
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
 * Free cuDSS linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.

 * parameters:
 *   context <-> pointer to cuDSS linear solver info
 *               (actual type: cs_sles_cudss_t  *)
 *----------------------------------------------------------------------------*/

void
cs_sles_cudss_free(void  *context);

/*----------------------------------------------------------------------------
 * Log sparse linear equation solver info.
 *
 * parameters:
 *   context  <-> pointer to cuDSS linear solver info
 *                (actual type: cs_sles_cudss_t  *)
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

void
cs_sles_cudss_log(const void  *context,
                  cs_log_t     log_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on cuDSS library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_library_info(cs_log_t  log_type);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set MPI communicator for cuDSS solver.
 *
 * The system is solved only on ranks with a non-null communicator or
 * if the associated communicator has less than 2 ranks.
 *
 * \param[in, out]  context  pointer to solver info and context
 * \param[in]       comm     MPI communicator for solving
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_set_mpi_comm(cs_sles_cudss_t  *context,
                           MPI_Comm          comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_SLES_CUDSS_H */
