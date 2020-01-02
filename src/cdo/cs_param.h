#ifndef __CS_PARAM_H__
#define __CS_PARAM_H__

/*============================================================================
 * Manage the definition/setting of a computation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 * \brief  Generic function pointer for defining a quantity at known locations
 *         (cells, faces, edges or vertices) with a function.  elt_ids is
 *         optional. If not NULL, the function works on a sub-list of
 *         elements. Moreover, it enables to fill retval with an indirection if
 *         compact is set to false
 *
 * \param[in]      n_elts   number of elements to consider
 * \param[in]      elt_ids  list of elements ids
 * \param[in]      compact  true:no indirection, false:indirection for retval
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] retval   result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_dof_func_t) (cs_lnum_t            n_elts,
                 const cs_lnum_t     *elt_ids,
                 bool                 compact,
                 void                *input,
                 cs_real_t           *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function which defines the evolution of a quantity according to the
 *         number of iteration already done, the current time and any structure
 *         given as a parameter
 *
 * \param[in]   time_iter   current number of iterations
 * \param[in]   time        value of the time at the end of the last iteration
 * \param[in]   input       pointer to a structure cast on-the-fly
 * \param[in]   retval      result of the evaluation
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_time_func_t) (int           time_iter,
                  double        time,
                  void         *input,
                  cs_real_t    *retval);

/* ================
 * ENUM definitions
 * ================ */

/*! \enum cs_param_space_scheme_t
 *  \brief Type of numerical scheme for the discretization in space
 *
 * \var CS_SPACE_SCHEME_LEGACY
 * Cell-centered Finite Volume Two-Point Flux
 *
 * \var CS_SPACE_SCHEME_CDOVB
 * CDO scheme with vertex-based positioning
 *
 * \var CS_SPACE_SCHEME_CDOVCB
 * CDO scheme with vertex+cell-based positioning
 *
 * \var CS_SPACE_SCHEME_CDOEB
 * CDO scheme with edge-based positioning
 *
 * \var CS_SPACE_SCHEME_CDOFB
 * CDO scheme with face-based positioning
 *
 * \var CS_SPACE_SCHEME_HHO_P0
 * Hybrid High Order (HHO) schemes
 * HHO scheme with face-based positioning (lowest order)
 *
 * \var CS_SPACE_SCHEME_HHO_P1
 * Hybrid High Order (HHO) schemes
 * HHO scheme with face-based positioning (k=1 up to order 3)
 *
 * \var CS_SPACE_SCHEME_HHO_P2
 * Hybrid High Order (HHO) schemes
 * HHO scheme with face-based positioning (k=2 up to order 4)
 */

typedef enum {

  CS_SPACE_SCHEME_LEGACY,
  CS_SPACE_SCHEME_CDOVB,
  CS_SPACE_SCHEME_CDOVCB,
  CS_SPACE_SCHEME_CDOEB,
  CS_SPACE_SCHEME_CDOFB,
  CS_SPACE_SCHEME_HHO_P0,
  CS_SPACE_SCHEME_HHO_P1,
  CS_SPACE_SCHEME_HHO_P2,

  CS_SPACE_N_SCHEMES

} cs_param_space_scheme_t;

/*! \enum cs_param_space_scheme_t
 *  \brief How is defined the degree of freedom
 *
 * \var CS_PARAM_REDUCTION_DERHAM
 *  Evaluation at vertices for potentials, integral along a line for
 *  circulations, integral across the normal component of a face for fluxes and
 *  integral inside a cell for densities
 *
 * \var CS_PARAM_REDUCTION_AVERAGE
 * Degrees of freedom are defined as the average on the element
 */

typedef enum {

  CS_PARAM_REDUCTION_DERHAM,
  CS_PARAM_REDUCTION_AVERAGE,

  CS_PARAM_N_REDUCTIONS

} cs_param_dof_reduction_t;

/*!
 * @name Numerical settings for time discretization
 * @{
 *
 * \enum cs_param_time_scheme_t
 *   Type of numerical scheme for the discretization in time
 *
 * \var CS_TIME_SCHEME_STEADY
 * No time scheme. Steady-state computation.
 *
 * \var CS_TIME_SCHEME_EULER_IMPLICIT
 * fully implicit (forward Euler/theta-scheme = 1)
 *
 * \var CS_TIME_SCHEME_EULER_EXPLICIT
 * fully explicit (backward Euler/theta-scheme = 0)
 *
 * \var CS_TIME_SCHEME_CRANKNICO
 * Crank-Nicolson (theta-scheme = 0.5)
 *
 * \var CS_TIME_SCHEME_THETA
 * theta-scheme
 */

typedef enum {

  CS_TIME_SCHEME_STEADY,
  CS_TIME_SCHEME_EULER_IMPLICIT,
  CS_TIME_SCHEME_EULER_EXPLICIT,
  CS_TIME_SCHEME_CRANKNICO,
  CS_TIME_SCHEME_THETA,

  CS_TIME_N_SCHEMES

} cs_param_time_scheme_t;

/*!
 * @}
 * @name Settings for the advection term
 * @{
 *
 * \enum cs_param_advection_form_t
 * Type of formulation for the advection term
 *
 * \var CS_PARAM_ADVECTION_FORM_CONSERV
 * conservative formulation (i.e. divergence of a flux)
 *
 * \var CS_PARAM_ADVECTION_FORM_NONCONS
 * non-conservative formation (i.e advection field times
 * the gradient)
 *
 * \var CS_PARAM_ADVECTION_FORM_SKEWSYM
 * skew-symmetric form
 *
 */

typedef enum {

  CS_PARAM_ADVECTION_FORM_CONSERV,
  CS_PARAM_ADVECTION_FORM_NONCONS,
  CS_PARAM_ADVECTION_FORM_SKEWSYM,

  CS_PARAM_N_ADVECTION_FORMULATIONS

} cs_param_advection_form_t;

/*! \enum cs_param_advection_scheme_t
 * Type of numerical scheme for the advection term
 *
 * \var CS_PARAM_ADVECTION_SCHEME_CENTERED
 * centered discretization
 *
 * \var CS_PARAM_ADVECTION_SCHEME_CIP
 * Continuous Interior Penalty discretization. Only available for
 * \ref CS_SPACE_SCHEME_CDOVCB
 *
 * \var CS_PARAM_ADVECTION_SCHEME_CIP_CW
 * Continuous Interior Penalty discretization. Only available for
 * \ref CS_SPACE_SCHEME_CDOVCB
 * Variant with a cellwise constant approximation of the advection field
 *
 * \var CS_PARAM_ADVECTION_SCHEME_MIX_CENTERED_UPWIND
 * Centered discretization with a portion between [0,1] of upwinding.
 * The portion is specified thanks to \ref CS_EQKEY_ADV_UPWIND_PORTION
 * If the portion is equal to 0, then one recovers
 * \ref CS_PARAM_ADVECTION_SCHEME_CENTERED. If the portion is equal to 1, then
 * one recovers CS_PARAM_ADVECTION_SCHEME_UPWIND
 *
 * \var CS_PARAM_ADVECTION_SCHEME_SAMARSKII
 * Weighting between an upwind and a centered discretization relying on the
 * Peclet number. Weighting function = Samarskii
 *
 * \var CS_PARAM_ADVECTION_SCHEME_SG
 * Weighting between an upwind and a centered discretization relying on the
 * Peclet number. Weighting function = Scharfetter-Gummel
 *
 * \var CS_PARAM_ADVECTION_SCHEME_UPWIND
 * Low order upwind discretization
 */

typedef enum {

  CS_PARAM_ADVECTION_SCHEME_CENTERED,
  CS_PARAM_ADVECTION_SCHEME_CIP,
  CS_PARAM_ADVECTION_SCHEME_CIP_CW,
  CS_PARAM_ADVECTION_SCHEME_MIX_CENTERED_UPWIND,
  CS_PARAM_ADVECTION_SCHEME_SAMARSKII,
  CS_PARAM_ADVECTION_SCHEME_SG,
  CS_PARAM_ADVECTION_SCHEME_UPWIND,

  CS_PARAM_N_ADVECTION_SCHEMES

} cs_param_advection_scheme_t;

/*!
 * @}
 * @name Settings for the boundary conditions
 * @{
 *
 * \enum cs_param_bc_type_t
 * Type of boundary condition to enforce.
 *
 * \var CS_PARAM_BC_HMG_DIRICHLET
 * Homogeneous Dirichlet boundary conditions. The value of a variable is set
 * to zero.
 *
 * \var CS_PARAM_BC_DIRICHLET
 * Dirichlet boundary conditions. The value of the variable is set to the user
 * requirements.
 *
 * \var CS_PARAM_BC_HMG_NEUMANN
 * Homogeneous Neumann conditions. The value of the flux of variable is set
 * to zero.
 *
 * \var CS_PARAM_BC_NEUMANN
 * Neumann conditions. The value of the flux of variable is set to the user
 * requirements.
 *
 * \var CS_PARAM_BC_ROBIN
 * Robin conditions.
 *
 * \var CS_PARAM_BC_SLIDING
 * Sliding conditions. Homogeneous Dirichlet for the normal component and
 * homogeneous Neumann for the tangential components. Only available for
 * vector-valued equations.
 *
 * \var CS_PARAM_BC_CIRCULATION
 * Set the tangential part of a vector-valued field. This is the part lying on
 * the boundary of a part of the computatinal domain. Nothing is prescribed for
 * the normal part of the vector-valued field.
 */

typedef enum {

  CS_PARAM_BC_HMG_DIRICHLET,
  CS_PARAM_BC_DIRICHLET,
  CS_PARAM_BC_HMG_NEUMANN,
  CS_PARAM_BC_NEUMANN,
  CS_PARAM_BC_ROBIN,
  CS_PARAM_BC_SLIDING,
  CS_PARAM_BC_CIRCULATION,
  CS_PARAM_N_BC_TYPES

} cs_param_bc_type_t;

/*! \enum cs_param_bc_enforce_t
 * Type of method for enforcing the boundary conditions. According to the type
 * of numerical scheme, some enforcements are not available.
 *
 * \var CS_PARAM_BC_ENFORCE_ALGEBRAIC
 * Weak enforcement of the boundary conditions (i.e. one keeps the degrees of
 * freedom in the algebraic system) with an algebraic manipulation of the
 * linear system
 *
 * \var CS_PARAM_BC_ENFORCE_PENALIZED
 * Weak enforcement of the boundary conditions (i.e. one keeps the degrees of
 * freedom in the algebraic system) with a penalization technique using a huge
 * value.
 *
 * \var CS_PARAM_BC_ENFORCE_WEAK_NITSCHE
 * Weak enforcement of the boundary conditions (i.e. one keeps the degrees of
 * freedom in the algebraic system) with a Nitsche-like penalization technique.
 * This technique does not increase the conditioning number as much as
 * \ref CS_PARAM_BC_ENFORCE_PENALIZED but is more computationally intensive.
 * The computation of the diffusion term should be activated with this choice.
 *
 * \var CS_PARAM_BC_ENFORCE_WEAK_SYM
 * Weak enforcement of the boundary conditions (i.e. one keeps the degrees of
 * freedom in the algebraic system) with a Nitsche-like penalization technique.
 * This variant enables to keep the symmetry of the algebraic contrary to
 * \ref CS_PARAM_BC_ENFORCE_WEAK_NITSCHE.
 * The computation of the diffusion term should be activated with this choice.
 */

typedef enum {

  CS_PARAM_BC_ENFORCE_ALGEBRAIC,
  CS_PARAM_BC_ENFORCE_PENALIZED,
  CS_PARAM_BC_ENFORCE_WEAK_NITSCHE,
  CS_PARAM_BC_ENFORCE_WEAK_SYM,

  CS_PARAM_N_BC_ENFORCEMENTS

} cs_param_bc_enforce_t;

/*!
 * @}
 * @name Settings for the linear solvers or SLES (Sparse Linear Equation Solver)
 * @{
 *
 * \enum cs_param_sles_class_t
 * \brief Class of iterative solvers to consider for solver the linear system
 *
 * \var CS_PARAM_SLES_CLASS_CS
 * Iterative solvers available in Code_Saturne
 *
 * \var CS_PARAM_SLES_CLASS_PETSC
 * Solvers available in PETSc
 *
 * \var CS_PARAM_SLES_CLASS_HYPRE
 * Solvers available in HYPRE through the PETSc library
 *
 * \var CS_PARAM_SLES_N_CLASSES
 */

typedef enum {

  CS_PARAM_SLES_CLASS_CS,
  CS_PARAM_SLES_CLASS_PETSC,
  CS_PARAM_SLES_CLASS_HYPRE,
  CS_PARAM_SLES_N_CLASSES

} cs_param_sles_class_t;

/*!
 * \enum cs_param_amg_type_t
 * Type of AMG (Algebraic MultiGrid) algorithm to use (either as a
 * preconditioner with or a solver).
 */

typedef enum {

  CS_PARAM_AMG_NONE,            /*!< No specified algorithm */
  CS_PARAM_AMG_HYPRE_BOOMER,    /*!< Boomer algorithm from Hypre library */
  CS_PARAM_AMG_PETSC_GAMG  ,    /*!< GAMG algorithm from PETSc */
  CS_PARAM_AMG_PETSC_PCMG  ,    /*!< preconditioned MG algorithm from PETSc */
  CS_PARAM_AMG_HOUSE_V,         /*!< In-house algorithm with V-cycle */
  CS_PARAM_AMG_HOUSE_K,         /*!< In-house algorithm with K-cycle */
  CS_PARAM_N_AMG_TYPES

} cs_param_amg_type_t;

/*!
 * \enum cs_param_precond_type_t
 * Type of preconditioner to use with the iterative solver.
 * Some of the mentionned preconditioners are available only if the PETSc
 * library is linked with code_saturne
 *
 * \var CS_PARAM_PRECOND_NONE
 * No preconditioner
 *
 * \var CS_PARAM_PRECOND_BJACOB
 * Block Jacobi
 *
 * \var CS_PARAM_PRECOND_AMG
 * Algebraic multigrid preconditioner (additional options may be set using
 * \ref \cs_param_amg_type_t)
 *
 * \var CS_PARAM_PRECOND_AMG_BLOCK
 * Algebraic multigrid preconditioner by block (useful in case of vector
 * valued variables)
 *
 * \var CS_PARAM_PRECOND_AS
 * Additive Schwarz preconditioner
 *
 * \var CS_PARAM_PRECOND_DIAG
 * Diagonal (also Jacobi) preconditioner. The cheapest one but not the most
 * efficient one.
 *
 * \var CS_PARAM_PRECOND_GKB_CG
 * Golub-Kahan Bidiagonalization solver used as a preconditioner. Only
 * useful if one has to solve a saddle-point system (such systems arise
 * when solving Stokes or Navier-Stokes in a fully couple manner).
 * Variant with CG as inner solver.
 *
 *
 * \var CS_PARAM_PRECOND_GKB_GMRES
 * Golub-Kahan Bidiagonalization solver used as a preconditioner. Only
 * useful if one has to solve a saddle-point system (such systems arise
 * when solving Stokes or Navier-Stokes in a fully couple manner).
 * Variant with GMRES as inner solver.
 *
 * \var CS_PARAM_PRECOND_ILU0
 * Incomplute LU factorization (fill-in coefficient set to 0)
 *
 * \var CS_PARAM_PRECOND_ICC0
 * Incomplute Cholesky factorization (fill-in coefficient set to 0). This is
 * variant of the ILU0 preconditioner dedicated to symmetric positive definite
 * system
 *
 * \var CS_PARAM_PRECOND_POLY1
 * Neumann polynomial preconditioning. Polynoms of order 1.
 *
 * \var CS_PARAM_PRECOND_POLY2
 * Neumann polynomial preconditioning. Polynoms of order 2.
 *
 * \var CS_PARAM_PRECOND_SSOR
 * Symmetric Successive OverRelaxations (can be seen as a symmetric
 * Gauss-Seidel preconditioner)
 */

typedef enum {

  CS_PARAM_PRECOND_NONE,

  CS_PARAM_PRECOND_BJACOB,
  CS_PARAM_PRECOND_AMG,
  CS_PARAM_PRECOND_AMG_BLOCK,
  CS_PARAM_PRECOND_AS,          /*!< Only with PETSc */
  CS_PARAM_PRECOND_DIAG,
  CS_PARAM_PRECOND_GKB_CG,
  CS_PARAM_PRECOND_GKB_GMRES,
  CS_PARAM_PRECOND_ILU0,        /*!< Only with PETSc */
  CS_PARAM_PRECOND_ICC0,        /*!< Only with PETSc*/
  CS_PARAM_PRECOND_POLY1,
  CS_PARAM_PRECOND_POLY2,
  CS_PARAM_PRECOND_SSOR,

  CS_PARAM_N_PRECOND_TYPES

} cs_param_precond_type_t;

/*!
 * \enum cs_param_itsol_type_t
 *  Type of solver to use to solve a linear system.
 *  Some of the mentionned solver are available only if the PETSc library is
 *  linked with code_saturne.
 *
 * \var CS_PARAM_ITSOL_NONE
 *  No iterative solver (equivalent to a "preonly" choice in PETSc)
 *
 * \var CS_PARAM_ITSOL_AMG
 *  Algebraic multigrid solver (additional options may be set using
 *  \ref \cs_param_amg_type_t)
 *
 * \var CS_PARAM_ITSOL_BICG
 *  Bi-Conjuguate gradient (useful for non-symmetric systems)
 *
 * \var CS_PARAM_ITSOL_BICGSTAB2
 *  Stabilized Bi-Conjuguate gradient (useful for non-symmetric systems)
 *
 * \var CS_PARAM_ITSOL_CG
 *  Conjuguate Gradient (solver of choice for symmetric positive definite
 *  systems)
 *
 * \var CS_PARAM_ITSOL_CR3
 * 3-layer conjugate residual (can handle non-symmetric systems)
 *
 * \var CS_PARAM_ITSOL_FCG
 * Flexible Conjuguate Gradient (variant of the CG when the preconditioner
 * may change from one iteration to another. For instance when using an AMG
 * as preconditioner)
 *
 * \var CS_PARAM_ITSOL_FGMRES
 * Flexible Generalized Minimal RESidual
 *
 * \var CS_PARAM_ITSOL_GAUSS_SEIDEL
 * Gauss-Seidel
 *
 * \var CS_PARAM_ITSOL_GKB_CG
 * Golub-Kahan Bidiagonalization algorithm. Useful for solving saddle-point
 * systems. The inner solver is a (flexible) CG solver.
 *
 * \var CS_PARAM_ITSOL_GKB_GMRES
 * Golub-Kahan Bidiagonalization algorithm. Useful for solving saddle-point
 * systems. The inner solver is a (flexible) GMRES solver.
 *
 * \var CS_PARAM_ITSOL_GMRES
 * Generalized Minimal RESidual
 *
 * \var CS_PARAM_ITSOL_JACOBI
 * Jacobi (diagonal relaxation)
 *
 * \var CS_PARAM_ITSOL_MINRES
 * Mininal residual algorithm
 *
 * \var CS_PARAM_ITSOL_MUMPS
 * MUMPS direct solver (LU factorization)
 *
 * \var CS_PARAM_ITSOL_MUMPS_LDLT
 * MUMPS direct solver (LDLT factorization also known as Cholesky factorization)
 *
 * \var CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL
 * Symmetric Gauss-Seidel
 *
 */

typedef enum {

  CS_PARAM_ITSOL_NONE,

  CS_PARAM_ITSOL_AMG,
  CS_PARAM_ITSOL_BICG,
  CS_PARAM_ITSOL_BICGSTAB2,
  CS_PARAM_ITSOL_CG,
  CS_PARAM_ITSOL_CR3,
  CS_PARAM_ITSOL_FCG,
  CS_PARAM_ITSOL_FGMRES,           /*!< Only with PETsc */
  CS_PARAM_ITSOL_GAUSS_SEIDEL,
  CS_PARAM_ITSOL_GKB_CG,
  CS_PARAM_ITSOL_GKB_GMRES,
  CS_PARAM_ITSOL_GMRES,            /*!< Only with PETsc */
  CS_PARAM_ITSOL_JACOBI,
  CS_PARAM_ITSOL_MINRES,           /*!< Only with PETsc */
  CS_PARAM_ITSOL_MUMPS,            /*!< Only with PETsc */
  CS_PARAM_ITSOL_MUMPS_LDLT,       /*!< Only with PETsc */
  CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL,
  CS_PARAM_N_ITSOL_TYPES

} cs_param_itsol_type_t;

/*!
 * \enum cs_param_resnorm_type_t
 * Way to renormalize (or not) the residual arising during the resolution of
 * linear systems
 */

typedef enum {

  CS_PARAM_RESNORM_NONE,           /*!< No renormalization  */
  CS_PARAM_RESNORM_NORM2_RHS,      /*!< Renormalization based on the Euclidean
                                        norm of the right-hand side */
  CS_PARAM_RESNORM_WEIGHTED_RHS,   /*!< Renormalization based on a weighted
                                        Euclidean norm of the right-hand side */
  CS_PARAM_RESNORM_FILTERED_RHS,   /*!< Renormalization based on an Euclidean
                                        norm of a selection of the right-hand
                                        side (penalized terms are filtered) */
  CS_PARAM_N_RESNORM_TYPES

} cs_param_resnorm_type_t;

/*!
 * \struct cs_param_sles_t
 * \brief Structure storing all metadata related to the resolution of a linear
 *        system with an iterative solver.
 */

typedef struct {

  bool                     setup_done;   /*!< SLES setup step has been done */
  int                      verbosity;    /*!< SLES verbosity */
  int                      field_id;     /*!< Field id related to a SLES
                                           By default, this is set to -1 */

  cs_param_sles_class_t    solver_class; /*!< class of SLES to consider  */
  cs_param_precond_type_t  precond;      /*!< type of preconditioner */
  cs_param_itsol_type_t    solver;       /*!< type of solver */
  cs_param_amg_type_t      amg_type;     /*!< type of AMG algorithm if needed */

  /*! \var resnorm_type
   *  normalized or not the norm of the residual used for the stopping criterion
   *  See \ref CS_EQKEY_ITSOL_RESNORM_TYPE for more details.
   */
  cs_param_resnorm_type_t  resnorm_type;
  int                      n_max_iter;   /*!< max. number of iterations */
  double                   eps;          /*!< stopping criterion on accuracy */

} cs_param_sles_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Return true if the space scheme has degrees of freedom on faces,
 *          otherwise false
 *
 * \param[in] scheme      type of space scheme
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_space_scheme_is_face_based(cs_param_space_scheme_t    scheme);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the space discretization scheme
 *
 * \param[in] scheme      type of space scheme
 *
 * \return the associated space scheme name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_space_scheme_name(cs_param_space_scheme_t    scheme);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the time discretization scheme
 *
 * \param[in] scheme      type of time scheme
 *
 * \return the associated time scheme name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_time_scheme_name(cs_param_time_scheme_t    scheme);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of boundary condition
 *
 * \param[in] type     type of boundary condition
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
 * \brief   Get the name of the type of algebraic multigrid (AMG)
 *
 * \param[in] type     type of AMG
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_amg_type_name(cs_param_amg_type_t   type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_H__ */
