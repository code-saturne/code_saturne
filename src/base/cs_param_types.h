#ifndef __CS_PARAM_TYPES_H__
#define __CS_PARAM_TYPES_H__

/*============================================================================
 * Set of definitions of structures and types used for setting a case
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @addtogroup Main DoF categories
 *
 * Main categories of Degrees of freedom (DoF) to consider for high-level
 * structures such as an interface, a range set, a matrix or an assembler
 *
 * Remark:
 * - scalar-valued HHO-P1 and vector-valued CDO-Fb shares the same structures
 * - face-based schemes or HHO 0th order schemes are equivalent
 *
 * @{
 */

#define CS_DOF_VTX_SCAL    0 /*!< Vertex or Vertex+Cell scalar-valued DoFs */
#define CS_DOF_VTX_VECT    1 /*!< Vertex or Vertex+Cell vector-valued DoFs */
#define CS_DOF_FACE_SCAL   2 /*!< Face or HHO 0th order scalar-valued DoFs */
#define CS_DOF_FACE_VECT   3 /*!< Face vector-valued DoFs */
#define CS_DOF_FACE_SCAP1  3 /*!< HHO 1st order polynomial scalar-valued DoFs */
#define CS_DOF_FACE_SCAP2  4 /*!< HHO 2nd order polynomial scalar-valued DoFs */
#define CS_DOF_FACE_VECP0  3 /*!< HHO 0th order polynomial vector-valued DoFs */
#define CS_DOF_FACE_VECP1  5 /*!< HHO 1st order polynomial vector-valued DoFs */
#define CS_DOF_FACE_VECP2  6 /*!< HHO 2nd order polynomial vector-valued DoFs */
#define CS_DOF_EDGE_SCAL   7 /*!< Edge-based scalar-valued DoFs */

#define CS_N_DOF_CASES     8

/*!
 *  @}
 *  @addtogroup scalar_params
 *
 * @{
 */

/*
 * Field property type
 */

/*! isotropic diffusion */

#define CS_ISOTROPIC_DIFFUSION (1 << 0)

/*! orthotropic diffusion */

#define CS_ORTHOTROPIC_DIFFUSION (1 << 1)

/*! diffusion by a left-multiplied symmetric 3x3 tensor */

#define CS_ANISOTROPIC_LEFT_DIFFUSION (1 << 2)

/*! diffusion by a right-multiplied symmetric 3x3 tensor */

#define CS_ANISOTROPIC_RIGHT_DIFFUSION (1 << 3)

/*! diffusion by a symmetric 3x3 tensor */

#define CS_ANISOTROPIC_DIFFUSION ((1 << 2) + (1 << 3))

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for an evaluation relying on an analytic
 *         function.
 *
 * For the calling function, elt_ids is optional. If not NULL, the coords
 * array should be accessed with an indirection. The same indirection can
 * be applied to fill retval if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids (in coords and retval)
 * \param[in]      coords        where ?
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_analytic_func_t) (cs_real_t            time,
                      cs_lnum_t            n_elts,
                      const cs_lnum_t     *elt_ids,
                      const cs_real_t     *coords,
                      bool                 dense_output,
                      void                *input,
                      cs_real_t           *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for computing a quantity at predefined
 *         locations such as degrees of freedom (DoF): cells, faces, edges or
 *         or vertices.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_dof_func_t) (cs_lnum_t            n_elts,
                 const cs_lnum_t     *elt_ids,
                 bool                 dense_output,
                 void                *input,
                 cs_real_t           *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function which defines the evolution of a quantity according to the
 *         current time and any structure given as a parameter
 *
 * \param[in]   time        value of the time at the end of the last iteration
 * \param[in]   input       NULL or pointer to a structure cast on-the-fly
 * \param[in]   retval      result of the evaluation
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_time_func_t) (double        time,
                  void         *input,
                  cs_real_t    *retval);

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
 * \var CS_SPACE_SCHEME_CDOCB
 * CDO scheme with cell-based positioning (potential at cells and fluxes at
 * faces)
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
 *
 * \var CS_SPACE_SCHEME_MACFB
 * MAC scheme with face-based positioning
 * faces)
 */

typedef enum {

  CS_SPACE_SCHEME_LEGACY,
  CS_SPACE_SCHEME_CDOVB,
  CS_SPACE_SCHEME_CDOVCB,
  CS_SPACE_SCHEME_CDOEB,
  CS_SPACE_SCHEME_CDOFB,
  CS_SPACE_SCHEME_CDOCB,
  CS_SPACE_SCHEME_HHO_P0,
  CS_SPACE_SCHEME_HHO_P1,
  CS_SPACE_SCHEME_HHO_P2,
  CS_SPACE_SCHEME_MACFB,

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
 *
 * \var CS_TIME_SCHEME_BDF2
 * theta-scheme
 */

typedef enum {

  CS_TIME_SCHEME_STEADY,
  CS_TIME_SCHEME_EULER_IMPLICIT,
  CS_TIME_SCHEME_EULER_EXPLICIT,
  CS_TIME_SCHEME_CRANKNICO,
  CS_TIME_SCHEME_THETA,
  CS_TIME_SCHEME_BDF2,

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
 * \ref CS_SPACE_SCHEME_CDOVCB space discretization.
 *
 * \var CS_PARAM_ADVECTION_SCHEME_CIP_CW
 * Continuous Interior Penalty discretization. Only available for
 * \ref CS_SPACE_SCHEME_CDOVCB space discretization.
 * Variant with a cellwise constant approximation of the advection field
 *
 * \var CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND
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
  CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND,
  CS_PARAM_ADVECTION_SCHEME_SAMARSKII,
  CS_PARAM_ADVECTION_SCHEME_SG,
  CS_PARAM_ADVECTION_SCHEME_UPWIND,

  CS_PARAM_N_ADVECTION_SCHEMES

} cs_param_advection_scheme_t;

/*! \enum cs_param_advection_strategy_t
 *  \brief Choice of how to handle the advection term in an equation
 *
 * \var CS_PARAM_ADVECTION_IMPLICIT_FULL
 * The advection term is implicitly treated. The non-linearity stemming from
 * this term has to be solved using a specific algorithm such as the Picard
 * (fixed-point) technique or more elaborated techniques.
 *
 * \var CS_PARAM_ADVECTION_IMPLICIT_LINEARIZED
 * The advection term is implicitly treated. The non-linearity stemming from
 * this term is simplified. Namely, one assumes a linearized advection. This is
 * equivalent to a one-step Picard technique.
 *
 * \var CS_PARAM_ADVECTION_EXPLICIT
 * The advection term is treated explicitly. One keeps the non-linearity
 * stemming from this term at the right hand-side. An extrapolation can be
 * used for the advection field (cf. \ref cs_param_advection_extrapol_t)
 */

typedef enum {

  CS_PARAM_ADVECTION_IMPLICIT_FULL,
  CS_PARAM_ADVECTION_IMPLICIT_LINEARIZED,
  CS_PARAM_ADVECTION_EXPLICIT,

  CS_PARAM_N_ADVECTION_STRATEGIES

} cs_param_advection_strategy_t;

/*! \enum cs_param_advection_extrapol_t
 *  \brief Choice of how to extrapolate the advection field in the advection
 *         term
 *
 * \var CS_PARAM_ADVECTION_EXTRAPOL_NONE
 * The advection field is not extrapolated. The last known advection field
 * is considered \f$ \phi^n \f$ if one computes \f$ u^(n+1) \f$ knowing
 * \f$ u^n. \f$ In case of a a non-linearity, this is \f$ \phi^(n+1,k) \f$
 * if one computes \f$ u^(n+1,k+1) \f$ knowing \f$  u^(n+1,k) \f$ and \f$ u^n.\f$
 * The initial step is such that \f$ u^(n+1,0) = u^n \f$.
 * This is a good choice with a first-order forward or backward time scheme.
 *
 * \var CS_PARAM_ADVECTION_EXTRAPOL_TAYLOR_2
 * The advection field is extrapolated with a 2nd order Taylor expansion
 * yielding \f$ \phi^extrap = 2 \phi^n - \phi^(n-1) \f$.
 * This corresponds to an estimation of \f$ \phi \f$ at n+1. Thus, this is
 * a good choice when associated with a BDF2 time scheme.
 *
 * \var CS_PARAM_ADVECTION_EXTRAPOL_ADAMS_BASHFORTH_2
 * The advection field is extrapolated with a 2nd order Adams-Bashforth
 * technique yielding \f$ \phi^extrap = 3/2 \phi^n - 1/2 \phi^(n-1) \f$.
 * This corresponds to an estimation of \f$ \phi \f$ at n+1/2. Thus, this
 *) is a good choice when associated with a Crank-Nilcolson time scheme.
 */

typedef enum {

  CS_PARAM_ADVECTION_EXTRAPOL_NONE,
  CS_PARAM_ADVECTION_EXTRAPOL_TAYLOR_2,
  CS_PARAM_ADVECTION_EXTRAPOL_ADAMS_BASHFORTH_2,

  CS_PARAM_N_ADVECTION_EXTRAPOLATIONS

} cs_param_advection_extrapol_t;

/*!
 * @}
 * @name Settings for the boundary conditions
 * @{
 *
 * \enum cs_param_bc_type_t
 * Type of boundary condition to enforce.
 *
 * \var CS_BC_HMG_DIRICHLET
 * Homogeneous Dirichlet boundary conditions. The value of a variable is set
 * to zero.
 *
 * \var CS_BC_DIRICHLET
 * Dirichlet boundary conditions. The value of the variable is set to the user
 * requirements.
 *
 * \var CS_BC_SYMMETRY
 * Homogeneous Neumann conditions. The value of the flux of variable is set
 * to zero.
 *
 * \var CS_BC_NEUMANN
 * Neumann conditions (along the face normal). The value of the flux of
 * variable is set to the user requirements.
 *
 * \var CS_BC_NEUMANN_FULL
 * Neumann conditions with full definition relative to all directions, i.e.
 * not only the normal direction. The value of the flux of variable is set to
 * the user requirements.
 *
 * \var CS_BC_ROBIN
 * Robin conditions.
 *
 * \var CS_BC_SYMMETRY
 * Sliding conditions. Homogeneous Dirichlet for the normal component and
 * homogeneous Neumann for the tangential components. Only available for
 * vector-valued equations.
 *
 * \var CS_BC_CIRCULATION
 * Set the tangential part of a vector-valued field. This is the part lying on
 * the boundary of a part of the computatinal domain. Nothing is prescribed for
 * the normal part of the vector-valued field.
 *
 * \var CS_BC_WALL_MODELLED
 * Hybrid way to enforce the value at a wall. Use a wall function to prescribe
 * the wall value knowing the value at the associated cell center for instance
 */

typedef enum {

  CS_BC_HMG_DIRICHLET = 0,
  CS_BC_DIRICHLET = 1,
  CS_BC_RADIATIVE_OUTLET = 2,
  CS_BC_NEUMANN = 3,
  CS_BC_SYMMETRY = 4,
  CS_BC_WALL_MODELLED = 5,
  CS_BC_ROUGH_WALL_MODELLED = 6,
  CS_BC_NEUMANN_FULL = 7,
  CS_BC_ROBIN = 9,
  CS_BC_CIRCULATION = 11,
  CS_BC_IMPOSED_TOT_FLUX = 13,
  CS_BC_GENERALIZED_SYM = 14,
  CS_BC_IMPOSED_EXCHANGE_COEF = 15,

  CS_PARAM_N_BC_TYPES

} cs_param_bc_type_t;

/* Compatibility macros */

#define  CS_PARAM_BC_HMG_DIRICHLET   CS_BC_HMG_DIRICHLET
#define  CS_PARAM_BC_DIRICHLET       CS_BC_DIRICHLET
#define  CS_PARAM_BC_HMG_NEUMANN     CS_BC_SYMMETRY
#define  CS_PARAM_BC_NEUMANN         CS_BC_NEUMANN
#define  CS_PARAM_BC_NEUMANN_FULL    CS_BC_NEUMANN_FULL
#define  CS_PARAM_BC_ROBIN           CS_BC_ROBIN
#define  CS_PARAM_BC_SLIDING         CS_BC_SYMMETRY
#define  CS_PARAM_BC_CIRCULATION     CS_BC_CIRCULATION
#define  CS_PARAM_BC_WALL_PRECRIBED  CS_BC_WALL_MODELLED

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

/*! \struct cs_param_convergence_t
 *  \brief Set of parameters to check the convergence (or the divergence) of an
 *         iterative process (tolerances or max. number of iterations)
 */

typedef struct {

/*!
 * \var atol
 * Absolute tolerance under which the iterative process is stopped
 *
 * \var rtol
 * Relative tolerance under which the iterative process is stopped
 *
 * \var dtol
 * Tolerance above which the iterative process is stated as "diverged".
 * Not used if < 0
 *
 * \var n_max_iter
 * Maximal number of iterations before stopping the iterative process
 */

  double               atol;
  double               rtol;
  double               dtol;
  int                  n_max_iter;

} cs_param_convergence_t;

/*!
 * @}
 * @name Settings for non-linear algorithms
 * @{
 *
 * \enum cs_param_nl_algo_t
 * \brief Class of non-linear iterative algorithm
 *
 * \var CS_PARAM_NL_ALGO_NONE
 * No specific algorithm. This mean a linear technique.
 *
 * \var CS_PARAM_NL_ALGO_PICARD
 * Picard (also called fixed point) algorithm
 *
 * \var CS_PARAM_NL_ALGO_ANDERSON
 * Anderson acceleration
 *
 * \var CS_PARAM_N_NL_ALGOS
 */

typedef enum {

  CS_PARAM_NL_ALGO_NONE,
  CS_PARAM_NL_ALGO_PICARD,
  CS_PARAM_NL_ALGO_ANDERSON,
  CS_PARAM_N_NL_ALGOS

} cs_param_nl_algo_t;

/*!
 * @}
 * @name Settings for the linear solvers or SLES (Sparse Linear Equation Solver)
 * @{
 */

/*!
 * \enum cs_param_solver_class_t
 * \brief Class of iterative solvers to consider for solver the linear system
 *
 * \var CS_PARAM_SOLVER_CLASS_CS
 * Iterative solvers available in code_saturne
 *
 * \var CS_PARAM_SOLVER_CLASS_HYPRE
 * Solvers available in HYPRE through the PETSc library
 *
 * \var CS_PARAM_SOLVER_CLASS_MUMPS
 * Solvers available with MUMPS (without the PETSc interface)
 *
 * \var CS_PARAM_SOLVER_CLASS_PETSC
 * Solvers available in PETSc. Please notice that
 * the MUMPS solver can be handled within PETSc if the installation of PETSc
 * includes the MUMPS library
 *
 * \var CS_PARAM_N_SOLVER_CLASSES
 */

typedef enum {

  CS_PARAM_SOLVER_CLASS_CS,
  CS_PARAM_SOLVER_CLASS_HYPRE,
  CS_PARAM_SOLVER_CLASS_MUMPS,
  CS_PARAM_SOLVER_CLASS_PETSC,

  CS_PARAM_N_SOLVER_CLASSES

} cs_param_solver_class_t;

/*!
 * \enum cs_param_precond_block_t
 * Type of preconditioning by block
 *
 * \var CS_PARAM_PRECOND_BLOCK_NONE
 * No block preconditioner is requested (default)
 *
 * \var CS_PARAM_PRECOND_BLOCK_DIAG
 * Only the diagonal blocks are considered in the preconditioner
 *
 * \var CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR
 * The diagonal blocks and the lower blocks are considered in the
 * preconditioner
 *
 * \var CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL
 * A symmetric Gauss-Seidel block preconditioning is considered
 * (cf. Y. Notay, "A new analysis of block preconditioners for saddle-point
 * problems" (2014), SIAM J. Matrix. Anal. Appl.)
 *
 * \var CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR
 * The diagonal blocks and the upper blocks are considered in the
 * preconditioner
 */

typedef enum {

  CS_PARAM_PRECOND_BLOCK_NONE,
  CS_PARAM_PRECOND_BLOCK_DIAG,
  CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR,
  CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL,
  CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR,

  CS_PARAM_N_PCD_BLOCK_TYPES

} cs_param_precond_block_t;

/*!
 * \enum cs_param_precond_type_t
 * Type of preconditioner to use with the iterative solver.
 * Some of the mentionned preconditioners are available only if the PETSc
 * library is linked with code_saturne
 *
 * \var CS_PARAM_PRECOND_NONE
 * No preconditioner
 *
 * \var CS_PARAM_PRECOND_BJACOB_ILU0
 * Block Jacobi with an ILU zero fill-in in each block
 *
 * \var CS_PARAM_PRECOND_BJACOB_SGS
 * Block Jacobi with a symmetric Gauss-Seidel in each block (rely on
 * Eisenstat's trick with PETSc)
 *
 * \var CS_PARAM_PRECOND_AMG
 * Algebraic multigrid preconditioner (additional options may be set using
 * \ref cs_param_amg_type_t)
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
 * \var CS_PARAM_PRECOND_GKB_GMRES
 * Golub-Kahan Bidiagonalization solver used as a preconditioner. Only
 * useful if one has to solve a saddle-point system (such systems arise
 * when solving Stokes or Navier-Stokes in a fully couple manner).
 * Variant with GMRES as inner solver.
 *
 * \var CS_PARAM_PRECOND_LU
 * LU factorization (direct solver) with PETSc
 *
 * \var CS_PARAM_PRECOND_ILU0
 * Incomplete LU factorization (fill-in coefficient set to 0)
 *
 * \var CS_PARAM_PRECOND_ICC0
 * Incomplete Cholesky factorization (fill-in coefficient set to 0). This is
 * variant of the ILU0 preconditioner dedicated to symmetric positive definite
 * system
 *
 * \var CS_PARAM_PRECOND_MUMPS
 * Sparse direct solver available with the MUMPS library and used as a
 * preconditioner
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
 *
 * \var CS_PARAM_PRECOND_HPDDM
 * Use HPDDM and Geneo approach for coarse space with PETSc
 */

typedef enum {

  CS_PARAM_PRECOND_NONE,

  CS_PARAM_PRECOND_AMG,
  CS_PARAM_PRECOND_BJACOB_ILU0, /*!< Only with PETSc */
  CS_PARAM_PRECOND_BJACOB_SGS,  /*!< Only with PETSc */
  CS_PARAM_PRECOND_DIAG,
  CS_PARAM_PRECOND_GKB_CG,    /*!< Only with PETSc */
  CS_PARAM_PRECOND_GKB_GMRES, /*!< Only with PETSc */
  CS_PARAM_PRECOND_LU,        /*!< Only with PETSc */
  CS_PARAM_PRECOND_ILU0,      /*!< Only with PETSc */
  CS_PARAM_PRECOND_ICC0,      /*!< Only with PETSc */
  CS_PARAM_PRECOND_MUMPS,     /*!< Only with MUMPS */
  CS_PARAM_PRECOND_HPDDM,     /*!< Only with PETSc */
  CS_PARAM_PRECOND_POLY1,
  CS_PARAM_PRECOND_POLY2,
  CS_PARAM_PRECOND_SSOR,

  CS_PARAM_N_PRECOND_TYPES

} cs_param_precond_type_t;

/*!
 * \enum cs_param_solver_type_t
 *  Type of solver used to solve a linear system.
 *  Some of the mentionned solver are available only if the PETSc library is
 *  linked with code_saturne.
 *
 * \var CS_PARAM_SOLVER_NONE
 *  No iterative solver (equivalent to a "preonly" choice in PETSc)
 *
 * \var CS_PARAM_SOLVER_AMG
 *  Algebraic multigrid solver (additional options may be set using
 *  \ref cs_param_amg_type_t)
 *
 * \var CS_PARAM_SOLVER_BICGS
 *  Bi-Conjugate gradient (useful for non-symmetric systems)
 *
 * \var CS_PARAM_SOLVER_BICGS2
 *  Stabilized Bi-Conjugate gradient (useful for non-symmetric systems)
 *
 * \var CS_PARAM_SOLVER_CG
 *  Conjugate Gradient (solver of choice for symmetric positive definite
 *  systems)
 *
 * \var CS_PARAM_SOLVER_CR3
 * 3-layer conjugate residual (can handle non-symmetric systems)
 *
 * \var CS_PARAM_SOLVER_FCG
 * Flexible Conjugate Gradient (variant of the CG when the preconditioner
 * may change from one iteration to another. For instance when using an AMG
 * as preconditioner)
 *
 * \var CS_PARAM_SOLVER_FGMRES
 * Flexible Generalized Minimal RESidual
 *
 * \var CS_PARAM_SOLVER_GAUSS_SEIDEL
 * Gauss-Seidel
 *
 * \var CS_PARAM_SOLVER_GCR
 * Generalized conjugate residual (flexible iterative solver for symmetric or
 * non-symmetric system)
 *
 * \var CS_PARAM_SOLVER_GMRES
 * Generalized Minimal RESidual
 *
 * \var CS_PARAM_SOLVER_JACOBI
 * Jacobi (diagonal relaxation)
 *
 * \var CS_PARAM_SOLVER_MUMPS
 * MUMPS sparse direct solver
 *
 * \var CS_PARAM_SOLVER_SYM_GAUSS_SEIDEL
 * Symmetric Gauss-Seidel
 *
 * \var CS_PARAM_SOLVER_USER_DEFINED
 * User-defined iterative solver. It relies on the implementation of the
 * the function cs_user_sles_it_solver()
 */

typedef enum {

  CS_PARAM_SOLVER_NONE,

  CS_PARAM_SOLVER_AMG,
  CS_PARAM_SOLVER_BICGS,
  CS_PARAM_SOLVER_BICGS2,
  CS_PARAM_SOLVER_CG,
  CS_PARAM_SOLVER_CR3,
  CS_PARAM_SOLVER_FCG,
  CS_PARAM_SOLVER_FGMRES,           /*!< Only with PETSc */
  CS_PARAM_SOLVER_GAUSS_SEIDEL,
  CS_PARAM_SOLVER_GCR,
  CS_PARAM_SOLVER_GMRES,            /*!< Only with PETSc */
  CS_PARAM_SOLVER_JACOBI,
  CS_PARAM_SOLVER_MUMPS,            /*!< Only with PETSc or MUMPS */
  CS_PARAM_SOLVER_SYM_GAUSS_SEIDEL,
  CS_PARAM_SOLVER_USER_DEFINED,

  CS_PARAM_N_SOLVER_TYPES

} cs_param_solver_type_t;

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
 * \enum cs_param_dotprod_type_t
 * Way to compute a dot product or its related norm
 */

typedef enum {

  CS_PARAM_DOTPROD_EUCLIDEAN,   /*!< Usual Eucliean 2-norm (no weight) */
  CS_PARAM_DOTPROD_CDO,         /*!< Weighted 2-norm from CDO quantities */

  CS_PARAM_N_DOTPROD_TYPES

} cs_param_dotprod_type_t;

/*!
 *  @}
 */

/*============================================================================
 * Global variables
 *============================================================================*/

/* Separation lines: header1, header2 (compatible with markdown), other */

extern const char cs_sep_h1[80];
extern const char cs_sep_h2[80];
extern const char cs_sepline[80];
extern const char cs_med_sepline[50];

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
 * \brief   Get the label associated to the advection formulation
 *
 * \param[in] adv_form      type of advection formulation
 *
 * \return the associated label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_form_name(cs_param_advection_form_t    adv_form);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label of the advection scheme
 *
 * \param[in] scheme      type of advection scheme
 *
 * \return the associated advection scheme label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_scheme_name(cs_param_advection_scheme_t    scheme);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label associated to the advection strategy
 *
 * \param[in] adv_stra      type of advection strategy
 *
 * \return the associated label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_strategy_name(cs_param_advection_strategy_t  adv_stra);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the label associated to the extrapolation used for the advection
 *         field
 *
 * \param[in] extrapol   type of extrapolation for the advection field
 *
 * \return the associated label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_extrapol_name(cs_param_advection_extrapol_t   extrapol);

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
 * \brief   Get the name of the non-linear algorithm
 *
 * \param[in] algo     type of algorithm
 *
 * \return the associated algorithm name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_nl_algo_name(cs_param_nl_algo_t   algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label (short name) of the non-linear algorithm
 *
 * \param[in] algo     type of algorithm
 *
 * \return the associated algorithm label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_nl_algo_label(cs_param_nl_algo_t   algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of dot product to apply
 *
 * \param[in] dp_type     type of dot product
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_dotprod_type_name(cs_param_dotprod_type_t   dp_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the name of the solver
 *
 * \param[in] solver  type of iterative solver
 *
 * \return the associated solver name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_solver_name(cs_param_solver_type_t  solver);

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
 * \brief   Get the name of the type of block preconditioning
 *
 * \param[in] type     type of block preconditioning
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_precond_block_name(cs_param_precond_block_t   type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_TYPES_H__ */
