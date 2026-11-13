#ifndef CS_SLES_AMGX_H
#define CS_SLES_AMGX_H

/*============================================================================
 * Sparse Linear Equation Solvers using AmgX
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
#include "alge/cs_matrix.h"
#include "alge/cs_sles.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * AmgX wrapper option flags
 */

/*! Use AMGX_comm_from_maps_1_ring instead of AMGX_distribution
  for parallel matrices when its use seems possible */
#define CS_SLES_AMGX_PREFER_COMM_FROM_MAPS     (1 << 0)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* AmgX linear solver context (opaque) */

typedef struct _cs_sles_amgx_t  cs_sles_amgx_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Define and associate an AmgX linear system solver for a given field or
 * equation name.
 *
 * If this system did not previously exist, it is added to the list of "known"
 * systems. Otherwise, its definition is replaced by the one defined here.
 *
 * This is a utility function: if finer control is needed, see cs_sles_define()
 * and cs_sles_amgx_create().
 *
 * Note that this function returns a pointer directly to the AmgX solver
 * management structure. This may be used to set further options.  If needed,
 * cs_sles_find() may be used to obtain a pointer to the matching cs_sles_t
 * container.
 */
/*----------------------------------------------------------------------------*/

extern "C" cs_sles_amgx_t *
cs_sles_amgx_define(int          f_id,
                    const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * Create AmgX linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix, the
 * matrix type will be forced to MATSHELL ("shell") regardless of the option
 * used.
 */
/*----------------------------------------------------------------------------*/

extern "C" cs_sles_amgx_t *
cs_sles_amgx_create(void);

/*----------------------------------------------------------------------------*/
/*
 * Create AmgX linear system solver info and context based on existing info and
 * context.
 */
/*----------------------------------------------------------------------------*/

extern "C" void *
cs_sles_amgx_copy(const void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Destroy AmgX linear system solver info and context.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_destroy(void  **context);

/*----------------------------------------------------------------------------*/
/*
 * Return the solver configuration for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration strings syntax.
 */
/*----------------------------------------------------------------------------*/

extern "C" const char *
cs_sles_amgx_get_config(void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Define the configuration for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration strings syntax.
 *
 * If this function is not called, a default configuration will be used.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_set_config(void        *context,
                        const char  *config);

/*----------------------------------------------------------------------------*/
/*
 * Return the name of the solver configuration file for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration file syntax.
 */
/*----------------------------------------------------------------------------*/

extern "C" const char *
cs_sles_amgx_get_config_file(void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Set the solver configuration file for an AmgX solver.
 *
 * Check the AmgX docummentation for configuration file syntax.
 *
 * If this function is not called, a default configuration will be used.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_set_config_file(void        *context,
                             const char  *path);

/*----------------------------------------------------------------------------*/
/*
 * Define whether an AmgX solver should use the device or host
 *
 * By default, the device will be used, but by callingg this function with
 * "use_device = false", only the host will be used.
 */
/*----------------------------------------------------------------------------*/

extern "C" bool
cs_sles_amgx_get_use_device(void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Define whether an AmgX solver should use the device or host.
 *
 * By default, the device will be used, but by callingg this function with
 * "use_device = false", only the host will be used.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_set_use_device(void  *context,
                            bool   use_device);

/*----------------------------------------------------------------------------*/
/*
 * Define additional AmgX solver usage flags
 *
 * By default, the device will be used, but by calling this function with
 * "use_device = false", only the host will be used.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_set_flags(void  *context,
                       int    flags);

/*----------------------------------------------------------------------------*/
/*
 * Query additional AmgX solver usage flags.
 */
/*----------------------------------------------------------------------------*/

extern "C" int
cs_sles_amgx_get_flags(void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Setup AmgX linear equation solver.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_setup(void               *context,
                   const char         *name,
                   const cs_matrix_t  *a,
                   int                 verbosity);

/*----------------------------------------------------------------------------*/
/*
 * Call AmgX linear equation solver.
 *
 * The precision, r_norm, and n_iter parameters are ignored here. The matching
 * configuration options should be set earlier, using the
 * cs_sles_amgx_set_config function
 */
/*----------------------------------------------------------------------------*/

extern "C" cs_sles_convergence_state_t
cs_sles_amgx_solve(void                *context,
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

/*----------------------------------------------------------------------------*/
/*
 * Free AmgX linear equation solver setup context.
 *
 * This function frees resolution-related data, such as buffers and
 * preconditioning but does not free the whole context, as info used for
 * logging (especially performance data) is maintained.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_free(void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Log sparse linear equation solver info.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_log(const void  *context,
                 cs_log_t     log_type);

/*----------------------------------------------------------------------------*/
/*
 * Print information on AmgX library.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_amgx_library_info(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/

#endif /* CS_SLES_AMGX_H */
