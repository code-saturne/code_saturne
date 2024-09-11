#ifndef __CS_CDOFB_MONOLITHIC_SLES_H__
#define __CS_CDOFB_MONOLITHIC_SLES_H__

/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
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
 * \brief Set pointers to shared structures
 *
 * \param[in] mesh     pointer to the mesh structure
 * \param[in] connect  pointer to additional CDO connectivities
 * \param[in] quant    pointer to additional CDO mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_sharing(const cs_mesh_t            *mesh,
                                      const cs_cdo_connect_t     *connect,
                                      const cs_cdo_quantities_t  *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the system helper for a CDO-Fb scheme solving the
 *        Navier-Stokes equation using a monolithic approach for the
 *        velocity-pressure coupling
 *
 * \param[in]      nsp      Navier-Stokes paremeters
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] sc       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_system_helper(const cs_navsto_param_t  *nsp,
                                            const cs_param_saddle_t  *saddlep,
                                            cs_cdofb_monolithic_t    *sc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the saddle solver and its context for a CDO-Fb scheme solving
 *        the Navier-Stokes equation using a monolithic approach for the
 *        velocity-pressure coupling
 *
 * \param[in]      nsp      set of parameters for the Navier-Stokes system
 * \param[in]      saddlep  parameters for solving a saddle-point problem
 * \param[in, out] sc       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_init_solver(const cs_navsto_param_t  *nsp,
                                     const cs_param_saddle_t  *saddlep,
                                     cs_cdofb_monolithic_t    *sc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Augmented Lagrangian-Uzawa algorithm.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_alu(const cs_navsto_param_t  *nsp,
                             cs_saddle_solver_t       *solver,
                             cs_real_t                *u_f,
                             cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach. The system is
 *        split into a velocity block and the (unassembled) divergence operator
 *        Block preconditioning using a Schur approximation on a Krylov solver
 *        such as the GCR or MINRES is available.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a saddle-point solver
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_block_krylov(const cs_navsto_param_t  *nsp,
                                      cs_saddle_solver_t       *solver,
                                      cs_real_t                *u_f,
                                      cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Golub-Kahan Bidiagonalization algorithm.
 *        In-house implementation. The PETSc implementation is also available
 *        but appears less efficient in our tests.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_gkb_inhouse(const cs_navsto_param_t  *nsp,
                                     cs_saddle_solver_t       *solver,
                                     cs_real_t                *u_f,
                                     cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation with a CDO face-based approach.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a saddle-point solver
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_full_system(const cs_navsto_param_t  *nsp,
                                     cs_saddle_solver_t       *solver,
                                     cs_real_t                *u_f,
                                     cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Notay's algebraic transformation.
 *        The full system is treated as one block and solved as it is.
 *        In this situation, PETSc or MUMPS are usually considered.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_notay(const cs_navsto_param_t  *nsp,
                               cs_saddle_solver_t       *solver,
                               cs_real_t                *u_f,
                               cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the Uzawa-CG algorithm.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_uzawa_cg(const cs_navsto_param_t  *nsp,
                                  cs_saddle_solver_t       *solver,
                                  cs_real_t                *u_f,
                                  cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a linear system arising from the discretization of the
 *        Navier-Stokes equation using a monolithic velocity-pressure coupling
 *        with a CDO face-based approach.
 *        Solve this system using the SIMPLE algorithm.
 *
 * \param[in]      nsp     set of parameters related to the Navier-Stokes eqs.
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] u_f     values of the velocity at faces (3 components)
 * \param[in, out] p_c     values of the pressure in cells
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_sles_simple(const cs_navsto_param_t  *nsp,
                                cs_saddle_solver_t       *solver,
                                cs_real_t                *u_f,
                                cs_real_t                *p_c);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_MONOLITHIC_SLES_H__ */
