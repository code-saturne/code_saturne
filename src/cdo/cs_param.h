#ifndef __CS_PARAM_H__
#define __CS_PARAM_H__

/*============================================================================
 * Manage the definition/setting of a computation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for an analytic function
 *         elt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *
 * \param[in]      time     when ?
 * \param[in]      n_elts   number of elements to consider
 * \param[in]      elt_ids  list of elements ids (to access coords and fill)
 * \param[in]      coords   where ?
 * \param[in]      compact  true:no indirection, false:indirection for filling
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] retval   result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_analytic_func_t) (cs_real_t            time,
                      cs_lnum_t            n_elts,
                      const cs_lnum_t     *elt_ids,
                      const cs_real_t     *coords,
                      bool                 compact,
                      void                *input,
                      cs_real_t           *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function which defines the time step according to the number of
 *         iteration already done, the current time and any structure given as
 *         a parameter
 *
 * \param[in]   time_iter   current number of iterations
 * \param[in]   time        value of the time at the end of the last iteration
 * \param[in]   input       pointer to a structure cast on-the-fly
 *
 * \return the value of the time step
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_timestep_func_t) (int       time_iter,
                      double    time,
                      void     *input);

/* ================
 * ENUM definitions
 * ================ */

/* Type of numerical scheme for the discretization in space */
typedef enum {

  CS_SPACE_SCHEME_LEGACY, /* Cell-centered Finite Volume Two-Point Flux */
  CS_SPACE_SCHEME_CDOVB,  /* CDO scheme with vertex-based positionning */
  CS_SPACE_SCHEME_CDOVCB, /* CDO scheme with vertex+cell-based positionning */
  CS_SPACE_SCHEME_CDOFB,  /* CDO cell-based scheme with hybridization */
  CS_SPACE_SCHEME_HHO_P0, /* Hybrid High Order scheme; P0 approx. of gradient */
  CS_SPACE_SCHEME_HHO_P1, /* Hybrid High Order scheme; P1 approx. of gradient */
  CS_SPACE_SCHEME_HHO_P2, /* Hybrid High Order scheme; P2 approx. of gradient */

  CS_SPACE_N_SCHEMES

} cs_param_space_scheme_t;

/* TIME SCHEME */
/* =========== */

/* Type of numerical scheme for the discretization in time */
typedef enum {

  /* fully implicit (forward Euler/theta-scheme = 1) */
  CS_TIME_SCHEME_IMPLICIT,
  /* fully explicit (backward Euler/theta-scheme = 0) */
  CS_TIME_SCHEME_EXPLICIT,
   /* Crank-Nicolson (theta-scheme = 0.5) */
  CS_TIME_SCHEME_CRANKNICO,
  /* theta-scheme */
  CS_TIME_SCHEME_THETA,
  /* Undefined */
  CS_TIME_N_SCHEMES

} cs_param_time_scheme_t;

/* ADVECTION OPERATOR PARAMETRIZATION */
/* ================================== */

typedef enum {

  CS_PARAM_ADVECTION_FORM_CONSERV,
  CS_PARAM_ADVECTION_FORM_NONCONS,
  CS_PARAM_N_ADVECTION_FORMULATIONS

} cs_param_advection_form_t;

/* Type of weight functions considered */
typedef enum {

  CS_PARAM_ADVECTION_SCHEME_CENTERED,
  CS_PARAM_ADVECTION_SCHEME_CIP,        // Only V+C CDO schemes
  CS_PARAM_ADVECTION_SCHEME_UPWIND,
  CS_PARAM_ADVECTION_SCHEME_SAMARSKII,
  CS_PARAM_ADVECTION_SCHEME_SG,
  CS_PARAM_N_ADVECTION_SCHEMES

} cs_param_advection_scheme_t;

/* BOUNDARY CONDITIONS */
/* =================== */

/* Physic-driven boundary */
typedef enum {

  CS_PARAM_BOUNDARY_WALL,
  CS_PARAM_BOUNDARY_INLET,
  CS_PARAM_BOUNDARY_OUTLET,
  CS_PARAM_BOUNDARY_SYMMETRY,
  CS_PARAM_N_BOUNDARY_TYPES

} cs_param_boundary_type_t;

/* Mathematical boundary conditions */
typedef enum {

  CS_PARAM_BC_HMG_DIRICHLET,
  CS_PARAM_BC_DIRICHLET,
  CS_PARAM_BC_HMG_NEUMANN,
  CS_PARAM_BC_NEUMANN,
  CS_PARAM_BC_ROBIN,
  CS_PARAM_N_BC_TYPES

} cs_param_bc_type_t;

/* Choice between different method to enforce the BCs */
typedef enum {

  CS_PARAM_BC_ENFORCE_STRONG,       /* Strong enforcement of the BCs */
  CS_PARAM_BC_ENFORCE_WEAK_PENA,    /* Weak enforcement with a strong
                                       penalization coefficient */
  CS_PARAM_BC_ENFORCE_WEAK_NITSCHE, /* Weak enforcement using Nitsche method */
  CS_PARAM_BC_ENFORCE_WEAK_SYM,     /* Weak enforcement using symmetric Nitsche
                                       method */
  CS_PARAM_N_BC_ENFORCEMENTS

} cs_param_bc_enforce_t;

/* ITERATIVE SOLVERS */
/* ================= */

typedef enum {

  CS_PARAM_PRECOND_NONE,    // No preconditioning
  CS_PARAM_PRECOND_DIAG,    // Diagonal preconditioning (also called Jacobi)
  CS_PARAM_PRECOND_BJACOB,  // Block Jacobi
  CS_PARAM_PRECOND_POLY1,   // Neumann polynomial preconditioning (Order 1)
  CS_PARAM_PRECOND_SSOR,    // Symmetric Successive OverRelaxations
  CS_PARAM_PRECOND_ILU0,    // Incomplete LU factorization
  CS_PARAM_PRECOND_ICC0,    // Incomplete Cholesky factorization
  CS_PARAM_PRECOND_AMG,     // Algebraic MultiGrid
  CS_PARAM_PRECOND_AS,      // Additive Schwarz method
  CS_PARAM_N_PRECOND_TYPES

} cs_param_precond_type_t;

/* Type of iterative solver to use to inverse the linear system */
typedef enum {

  CS_PARAM_ITSOL_JACOBI,    // Jacobi
  CS_PARAM_ITSOL_CG,        // Conjuguate Gradient
  CS_PARAM_ITSOL_BICG,      // Bi-Conjuguate gradient
  CS_PARAM_ITSOL_BICGSTAB2, // Stabilized Bi-Conjuguate gradient
  CS_PARAM_ITSOL_CR3,       // 3-layer conjugate residual
  CS_PARAM_ITSOL_GMRES,     // Generalized Minimal RESidual
  CS_PARAM_ITSOL_AMG,       // Algebraic MultiGrid
  CS_PARAM_N_ITSOL_TYPES

} cs_param_itsol_type_t;

/* Description of the algorithm used to solve an equation */
typedef struct {

  cs_param_precond_type_t  precond; /* type of preconditioner */
  cs_param_itsol_type_t    solver;  /* type of solver */

  int          n_max_iter;          /* max. number of iterations */
  double       eps;                 /* stopping criterion on accuracy */
  bool         resid_normalized;    /* normalized or not the norm of the
                                       residual used for the stopping
                                       criterion */
} cs_param_itsol_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the solver
 *
 * \param[in] solver     type of iterative solver
 *
 * \return the associated solver name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_solver_name(cs_param_itsol_type_t  solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the preconditioner
 *
 * \param[in] precond     type of preconditioner
 *
 * \return the associated preconditioner name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_precond_name(cs_param_precond_type_t  precond);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of boundary condition
 *
 * \param[in] bc          type of boundary condition
 *
 * \return the associated bc name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_name(cs_param_bc_type_t  bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of enforcement of the boundary condition
 *
 * \param[in] type          type of enforcement of boundary conditions
 *
 * \return the associated name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_enforcement_name(cs_param_bc_enforce_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the domain boundary condition
 *          This name is also used as a name for zone definition
 *
 * \param[in] type     type of boundary
 *
 * \return the associated boundary name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_boundary_domain_name(cs_param_boundary_type_t  type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_H__ */
