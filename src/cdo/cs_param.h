#ifndef __CS_PARAM_H__
#define __CS_PARAM_H__

/*============================================================================
 * Manage the definition/setting of a computation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_mesh.h"
#include "cs_cdo.h"
#include "cs_quadrature.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* DISCRETE HODGE OPERATORS */
/* ======================== */

typedef enum {

  CS_PARAM_HODGE_TYPE_VPCD, // from primal vertices to dual cells
  CS_PARAM_HODGE_TYPE_EPFD, // from primal edges to dual faces
  CS_PARAM_HODGE_TYPE_FPED, // from primal faces to dual edges
  CS_PARAM_HODGE_TYPE_EDFP, // from dual edges to primal faces
  CS_PARAM_HODGE_TYPE_CPVD, // from primal cells to dual vertices
  CS_PARAM_HODGE_TYPE_VC,   // primal vertices + primal cells
  CS_PARAM_N_HODGE_TYPES

} cs_param_hodge_type_t;

typedef enum {

  CS_PARAM_HODGE_ALGO_VORONOI, // Under othogonality condition gives a diag. op.
  CS_PARAM_HODGE_ALGO_WBS,     // WBS: Whitney Barycentric Subdivision
  CS_PARAM_HODGE_ALGO_COST,    // COST: COnsistency & STabilization splitting
  CS_PARAM_HODGE_ALGO_AUTO,    /* Switch between the previous algo. according to
                                  the type of cell and its property */
  CS_PARAM_N_HODGE_ALGOS

} cs_param_hodge_algo_t;

typedef struct {

  bool   is_unity; /* No associated property. Property is equalt to the unity */
  bool   is_iso;   /* Associated property is isotroopic ? */
  bool   inv_pty;  /* Definition based on the property or its inverse */

  cs_param_hodge_type_t   type;  /* type of discrete Hodge operator */
  cs_param_hodge_algo_t   algo;  /* type of algorithm used to build this op. */
  double                  coef;  /* Value of the stabilization parameter
                                    if the COST algo. is used, otherwise 0. */

} cs_param_hodge_t;

/* TIME SCHEME */
/* =========== */

/* Type of numerical scheme for the discretization in time */
typedef enum {

  CS_TIME_SCHEME_IMPLICIT,  // fully implicit (forward Euler/theta-scheme = 1)
  CS_TIME_SCHEME_EXPLICIT,  // fully explicit (backward Euler/theta-scheme = 0)
  CS_TIME_SCHEME_CRANKNICO, // Crank-Nicolson (theta-scheme = 0.5)
  CS_TIME_SCHEME_THETA,     // theta-scheme
  CS_TIME_N_SCHEMES

} cs_time_scheme_t;

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

/* Choice on the type of algo. used to set the argument used in
   the weight function */
typedef enum {

  CS_PARAM_ADVECTION_WEIGHT_FLUX,
  CS_PARAM_ADVECTION_WEIGHT_XEXC,
  CS_PARAM_N_ADVECTION_WEIGHTS

} cs_param_advection_weight_t;

/* Set of options for building a contraction operator (also called interior
   product) which is closely related to the advection operator */
typedef struct {

  cs_param_advection_form_t     formulation; // conservative or not
  cs_param_advection_scheme_t   scheme;
  cs_param_advection_weight_t   weight_criterion;

} cs_param_advection_t;

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

  cs_param_precond_type_t  precond; // type of preconditioner
  cs_param_itsol_type_t    solver;  // type of solver

  int          n_max_iter;          // max. number of iterations
  double       eps;                 // stopping criterion on accuracy

  int          output_freq;         // frequencdy of output into listing
  bool         resid_normalized;    /* normalized or not the norm of the
                                       residual used for the stopping criterion
                                    */
} cs_param_itsol_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of algorithm related to a discrete Hdoge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the algorithm
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_algo_name(const cs_param_hodge_t   h_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the type of discrete Hodge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_type_name(const cs_param_hodge_t   h_info);

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

END_C_DECLS

#endif /* __CS_PARAM_H__ */
