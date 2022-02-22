#ifndef __CS_SLES_MUMPS_H__
#define __CS_SLES_MUMPS_H__

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

/*----------------------------------------------------------------------------
 * MUMPS headers
 *----------------------------------------------------------------------------*/

#include <dmumps_c.h>
#include <smumps_c.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo_perio.h"
#include "cs_matrix.h"
#include "cs_param_sles.h"
#include "cs_time_plot.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_SLES_MUMPS_JOB_INIT -1

/* Set of macros defined in order to match MUMPS documentation (because of the
 * difference between C/FORTRAN programming language
 */
#define ICNTL(I)   icntl[(I)-1]
#define CNTL(I)    cntl[(I)-1]
#define INFOG(I)   infog[(I)-1]
#define INFO(I)    info[(I)-1]
#define RINFOG(I)  rinfog[(I)-1]
#define RINFO(I)   rinfo[(I)-1]

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for user settings of a MUMPS solver.
 *        This function is called at the end of the setup stage.
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * \param[in]      slesp     pointer to the related cs_param_sles_t structure
 * \param[in, out] context   pointer to optional (untyped) value or structure
 * \param[in, out] dmumps    pointer to DMUMPS_STRUC_C (double-precision)
 * \param[in, out] smumps    pointer to SMUMPS_STRUC_C (single-precision)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_sles_mumps_setup_hook_t) (const cs_param_sles_t   *slesp,
                              void                    *context,
                              DMUMPS_STRUC_C          *dmumps,
                              SMUMPS_STRUC_C          *smumps);

/* MUMPS solver context (opaque) */

typedef struct _cs_sles_mumps_t  cs_sles_mumps_t;

/*============================================================================
 *  Global variables
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
                        SMUMPS_STRUC_C          *smumps);

/*=============================================================================
 * Public function prototypes
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
                     void                          *context);

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
                     void                        *context);

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
cs_sles_mumps_copy(const void   *context);

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
cs_sles_mumps_free(void  *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy MUMPS linear system solver info and context.
 *
 * \param[in, out]  context  pointer to sparse direct solver info and context
 *                           (actual type: cs_sles_mumps_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_mumps_destroy(void   **context);

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
                    int                 verbosity);

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
                    void                *aux_vectors);

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
                  cs_log_t     log_type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_MUMPS_H__ */
