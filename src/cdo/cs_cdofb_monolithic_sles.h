#ifndef __CS_CDOFB_MONOLITHIC_SLES_H__
#define __CS_CDOFB_MONOLITHIC_SLES_H__

/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_monolithic_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdofb_monolithic_sles_t structure
 *
 * \param[in] n_faces     number of faces (interior + border)
 * \param[in] n_cells     number of cells
 *
 * \return a pointer to a newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdofb_monolithic_sles_t *
cs_cdofb_monolithic_sles_create(cs_lnum_t    n_faces,
                                cs_lnum_t    n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a part of the structure
 *
 * \param[in, out]  msles   pointer to the structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_clean(cs_cdofb_monolithic_sles_t   *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free memory related to cs_cdofb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_free(cs_cdofb_monolithic_sles_t   **p_msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to shared structures
 *
 * \param[in]  connect  pointer to cdo connectivities
 * \param[in]  quant    pointer to additional mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_sharing(const cs_cdo_connect_t        *connect,
                                      const cs_cdo_quantities_t     *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free if needed structure(s) associated CDO face-based schemes with
 *         a monolithic velocity-pressure coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage.
 *         nsp is not declared as const to avoid compilation warnings but
 *         it should be modified at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_set_sles(cs_navsto_param_t    *nsp,
                             void                 *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from the discretization of the
 *         Navier-Stokes equation with a CDO face-based approach.
 *         The full system is treated as one block and solved as it is.
 *         In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_solve(const cs_navsto_param_t       *nsp,
                          const cs_equation_param_t     *eqp,
                          const cs_cdo_system_helper_t  *sh,
                          cs_param_sles_t               *slesp,
                          cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach. The system is
 *        split into blocks to enable more efficient preconditioning
 *        techniques. The main iterative solver is a Krylov solver such as GCR,
 *        or MINRES
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_krylov_block_precond(const cs_navsto_param_t       *nsp,
                                         const cs_equation_param_t     *eqp,
                                         const cs_cdo_system_helper_t  *sh,
                                         cs_param_sles_t               *slesp,
                                         cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the GKB algorithm to solve the saddle-point problem arising
 *         from CDO-Fb schemes for Stokes and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_gkb_solve(const cs_navsto_param_t       *nsp,
                              const cs_equation_param_t     *eqp,
                              const cs_cdo_system_helper_t  *sh,
                              cs_param_sles_t               *slesp,
                              cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the preconditioned Uzawa-CG algorithm to solve the saddle-point
 *         problem arising from CDO-Fb schemes for Stokes, Oseen and
 *         Navier-Stokes with a monolithic coupling
 *         This algorithm is based on Koko's paper "Uzawa conjugate gradient
 *         method for the Stokes problem: Matlab implementation with P1-iso-P2/
 *         P1 finite element"
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_cg_solve(const cs_navsto_param_t       *nsp,
                                   const cs_equation_param_t     *eqp,
                                   const cs_cdo_system_helper_t  *sh,
                                   cs_param_sles_t               *slesp,
                                   cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the Uzawa algorithm with an Augmented Lagrangian (ALU) technique
 *         in an incremental way to solve the saddle-point problem arising from
 *         CDO-Fb schemes for Stokes, Oseen and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      sh       pointer to a cs_cdo_system_helper_t structure
 * \param[in, out] slesp    pointer to a set of parameters to drive the SLES
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_al_incr_solve(const cs_navsto_param_t       *nsp,
                                        const cs_equation_param_t     *eqp,
                                        const cs_cdo_system_helper_t  *sh,
                                        cs_param_sles_t               *slesp,
                                        cs_cdofb_monolithic_sles_t    *msles);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_MONOLITHIC_SLES_H__ */
