#ifndef __CS_CDOCB_SCALEQ_SLES_H__
#define __CS_CDOCB_SCALEQ_SLES_H__

/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO cell-based schemes with a scaleq equations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cdo/cs_cdocb_priv.h"

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
 * \brief Set pointers to shared structures
 *
 * \param[in] mesh     pointer to the mesh structure
 * \param[in] connect  pointer to additional CDO connectivities
 * \param[in] quant    pointer to additional CDO mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_sles_init_sharing(const cs_mesh_t            *mesh,
                                  const cs_cdo_connect_t     *connect,
                                  const cs_cdo_quantities_t  *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the system helper for a CDO-Cb scheme solving a scalar-valued
 *        equation (saddle-point system)
 *
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] eqc      pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_sles_init_system_helper(const cs_param_saddle_t  *saddlep,
                                        cs_cdocb_scaleq_t        *eqc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the saddle solver and its context for a CDO-Cb scheme solving
 *        a scalar-valued equation (saddle-point problem)
 *
 * \param[in]      eqp      set of equation parameters
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] eqc      pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_sles_init_solver(const cs_equation_param_t  *eqp,
                                 const cs_param_saddle_t    *saddlep,
                                 cs_cdocb_scaleq_t          *eqc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes. Solve
 *        the saddle-point system using the Augmented Lagrangian-Uzawa
 *        algorithm.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_alu(cs_saddle_solver_t  *solver,
                         cs_real_t           *flux,
                         cs_real_t           *pot);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes. The
 *        system is split into a flux block and the (unassembled) divergence
 *        operator.
 *        Block preconditioning using a Schur approximation on a Krylov solver
 *        such as the GCR or MINRES is available.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_block_krylov(cs_saddle_solver_t  *solver,
                                  cs_real_t           *flux,
                                  cs_real_t           *pot);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the saddle-point linear system arising from the discretization
 *        of the scalar-valued CDO cell-based scheme.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in, out] solver  pointer to a saddle-point solver
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_full_system(cs_saddle_solver_t  *solver,
                                 cs_real_t           *flux,
                                 cs_real_t           *pot);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes.
 *        Solve this system using the Golub-Kahan Bidiagonalization algorithm.
 *        In-house implementation. The PETSc implementation is also available
 *        but appears less efficient in our tests.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_gkb_inhouse(cs_saddle_solver_t  *solver,
                                 cs_real_t           *flux,
                                 cs_real_t           *pot);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes.
 *        Solve this system using the Notay's algebraic transformation.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_notay(cs_saddle_solver_t  *solver,
                           cs_real_t           *flux,
                           cs_real_t           *pot);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of a
 *        scalar-valed equation stemming from CDO cell-based schemes.
 *        Solve this system using the Uzawa-CG algorithm.
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] flux    values of the flux at faces (scalar)
 * \param[in, out] pot     values of the potential in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdocb_scaleq_sles_uzawa_cg(cs_saddle_solver_t  *solver,
                              cs_real_t           *flux,
                              cs_real_t           *pot);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOCB_SCALEQ_SLES_H__ */
