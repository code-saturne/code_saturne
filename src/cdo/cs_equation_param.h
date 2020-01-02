#ifndef __CS_EQUATION_PARAM_H__
#define __CS_EQUATION_PARAM_H__

/*============================================================================
 * Header to handle specific settings related to a cs_equation_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_property.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_equation_param.h

  \brief Structure and routines handling the specific settings related
         to a cs_equation_t structure

*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @name Flags specifying which term is needed for an equation.
 * @{
 *
 * \def CS_EQUATION_LOCKED
 * \brief Parameters for setting an equation are not modifiable anymore
 *
 * \def CS_EQUATION_UNSTEADY
 * \brief Unsteady term is needed
 *
 * \def CS_EQUATION_CONVECTION
 * \brief Convection term is needed
 *
 * \def CS_EQUATION_DIFFUSION
 * \brief Diffusion term is needed. A scalar-/vector-valued Laplacian
 * with div .grad
 *
 * \def CS_EQUATION_CURLCURL
 * \brief The term related to the curl-curl operator is needed
 *
 * \def CS_EQUATION_GRADDIV
 * \brief The term related to the grad-div operator is needed
 *
 * \def CS_EQUATION_REACTION
 * \brief Reaction term is needed
 *
 * \def CS_EQUATION_FORCE_VALUES
 * \brief Add an algebraic manipulation to set the value of a given set of
 *       interior degrees of freedom
 *
 */

#define CS_EQUATION_LOCKED        (1 <<  0)  /*   1 */
#define CS_EQUATION_UNSTEADY      (1 <<  1)  /*   2 */
#define CS_EQUATION_CONVECTION    (1 <<  2)  /*   4 */
#define CS_EQUATION_DIFFUSION     (1 <<  3)  /*   8 */
#define CS_EQUATION_CURLCURL      (1 <<  4)  /*  16 */
#define CS_EQUATION_GRADDIV       (1 <<  5)  /*  32 */
#define CS_EQUATION_REACTION      (1 <<  6)  /*  64 */
#define CS_EQUATION_FORCE_VALUES  (1 <<  7)  /* 128 */

/*!
 * @}
 * @name Flags specifying which extra operation is needed for an equation.
 * @{
 *
 * \def CS_EQUATION_POST_BALANCE
 * \brief Compute and postprocess the equation balance
 *
 * \def CS_EQUATION_POST_PECLET
 * \brief Compute and postprocess the Peclet number
 *
 * \def CS_EQUATION_POST_UPWIND_COEF
 * \brief Postprocess the value of the upwinding coefficient
 *
 * \def CS_EQUATION_POST_NORMAL_FLUX
 * \brief Postprocess the value of the normal flux across the boundary faces
 *
 */

#define CS_EQUATION_POST_BALANCE     (1 << 0) /* 1 */
#define CS_EQUATION_POST_PECLET      (1 << 1) /* 2 */
#define CS_EQUATION_POST_UPWIND_COEF (1 << 2) /* 4 */
#define CS_EQUATION_POST_NORMAL_FLUX (1 << 3) /* 8 */

/*!
 * @}
 * @name Flags to handle the enforcement of degrees of freedom (DoFs)
 * @{
 *
 * \def CS_EQUATION_ENFORCE_BY_CELLS
 * \brief Definition of a selection of DoFs to enforce using a cell selection
 *
 * \def CS_EQUATION_ENFORCE_BY_DOFs
 * \brief Definition of a selection of DoFs
 *
 * \def CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE
 * \brief Assign to all the selected DoFs the same value. This value is stored
 *        in enforcement_ref_value
 *
 */

#define CS_EQUATION_ENFORCE_BY_CELLS            (1 << 0) /* 1 */
#define CS_EQUATION_ENFORCE_BY_DOFS             (1 << 1) /* 2 */
#define CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE  (1 << 2) /* 4 */

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_equation_type_t
 *  \brief Type of equations managed by the solver
 *
 * \var CS_EQUATION_TYPE_USER
 * User-defined equation
 *
 * \var CS_EQUATION_TYPE_GROUNDWATER
 * Equation related to the groundwater flow module
 *
 * \var CS_EQUATION_TYPE_MAXWELL
 * Equation related to the Maxwell module
 *
 * \var CS_EQUATION_TYPE_NAVSTO
 * Equation related to the resolution of the Navier-Stokes system
 * - Example: momentum, prediction, correction, energy...
 *
 * \var CS_EQUATION_TYPE_THERMAL
 * Equation related to the heat transfer
 *
 * \var CS_EQUATION_TYPE_PREDEFINED
 * Predefined equation (most part of the setting is already done)
 * - Example: equation for the wall distance or ALE
 */

typedef enum {

  CS_EQUATION_TYPE_USER,
  CS_EQUATION_TYPE_GROUNDWATER,
  CS_EQUATION_TYPE_MAXWELL,
  CS_EQUATION_TYPE_NAVSTO,
  CS_EQUATION_TYPE_THERMAL,
  CS_EQUATION_TYPE_PREDEFINED,
  CS_EQUATION_N_TYPES

} cs_equation_type_t;

/*! \struct cs_equation_param_t
 *  \brief Set of parameters to handle an unsteady convection-diffusion-reaction
 *         equation with term sources
 */

typedef struct {

  /*!
   * @name General settings
   * @{
   */
  char *restrict       name;           /*!< name of the equation */
  cs_equation_type_t   type;           /*!< type of equation: predefined... */
  int                  dim;            /*!< Dimension of the unknown */
  int                  verbosity;      /*!< Level of detail for output */

  /*! \var flag
   *  Flag to know if unsteady or diffusion or convection or reaction
   *  or source term are activated or not
   */
  cs_flag_t                  flag;

  /*! \var process_flag
   *  Flag to known if predefined post-treatments such as Peclet,
   *  are requested
   */
  cs_flag_t                  process_flag;

  /* Numerical settings */
  cs_param_space_scheme_t    space_scheme;  /*!< Space discretization scheme */
  cs_param_dof_reduction_t   dof_reduction; /*!< How is defined DoF */

  /*! \var space_poly_degree
   * Maximum degree of the polynomial basis
   */
  int                        space_poly_degree;

  /*!
   * @}
   * @name Settings for the boundary conditions
   * @{
   *
   * \var default_bc
   * Default boundary condition related to this equation. Valid choices:
   * - \ref CS_PARAM_BC_HMG_NEUMANN
   * - \ref CS_PARAM_BC_HMG_DIRICHLET
   *
   * \var n_bc_defs
   * Number of boundary conditions which are defined for this equation
   *
   * \var bc_defs
   * Pointers to the definitions of the boundary conditions
   *
   * \var default_enforcement
   * Type of strategy to enforce an essential boundary conditions (Dirichlet for
   * instance) when no predefined strategy is prescribed.  See \ref
   * cs_param_bc_enforce_t for more details.
   *
   * \var strong_pena_bc_coeff
   * Value of the penalization coefficient used to enforce the Dirichlet
   * boundary conditions when \ref CS_PARAM_BC_ENFORCE_PENALIZED is set. This
   * value should be sufficiently large in order to neglect off-diagonal terms.
   *
   * \var weak_pena_bc_coeff
   * Value of the penalization coefficient used to enforce the Dirichlet
   * boundary condition when \ref CS_PARAM_BC_ENFORCE_WEAK_NITSCHE or \ref
   * CS_PARAM_BC_ENFORCE_WEAK_SYM is set. This two latter strategies have a
   * lesser influence on the conditioning number of the linear system than the
   * choice \ref CS_PARAM_BC_ENFORCE_PENALIZED
   *
   */

  cs_param_bc_type_t            default_bc;
  int                           n_bc_defs;
  cs_xdef_t                   **bc_defs;

  cs_param_bc_enforce_t         default_enforcement;
  cs_real_t                     strong_pena_bc_coeff;
  cs_real_t                     weak_pena_bc_coeff;

  /*!
   * @}
   * @name Numerical settings for the time-dependent parameters
   * @{
   *
   * \var n_ic_defs
   * Number of definitions for setting the intial condition
   *
   * \var ic_defs
   * List of pointers to the definition of the inititial condition
   */

  int                           n_ic_defs;
  cs_xdef_t                   **ic_defs;

  /*!
   * @}
   * @name Numerical settings for the time-dependent parameters
   * @{
   *
   * \var do_lumping
   * Set to true or false. Activate several actions:
   * Perform a mass lumping of the matrices related to the time and reaction
   * discretization. All source terms are evaluated using a barycentric
   * quadrature.
   */

  bool                          do_lumping;

  /*!
   * \var time_hodge
   * Set of parameters for the discrete Hodge operator related to the unsteady
   * term
   *
   * \var time_property
   * Pointer to the \ref cs_property_t structure related to the unsteady term
   *
   * \var time_scheme
   * Type of numerical scheme used for the time discretization
   *
   * \var theta
   * Value of the coefficient for a theta scheme (between 0 and 1)
   *
   */

  cs_param_hodge_t              time_hodge;
  cs_property_t                *time_property;
  cs_param_time_scheme_t        time_scheme;
  cs_real_t                     theta;

  /*!
   * @}
   * @name Numerical settings for the diffusion (Laplacian div-grad) term
   * @{
   *
   * \var diffusion_hodge
   * Set of parameters for the discrete Hodge operator related to the diffusion
   *
   * \var diffusion_property
   * Pointer to the property related to the diffusion term
   */

  cs_param_hodge_t              diffusion_hodge;
  cs_property_t                *diffusion_property;

  /*!
   * @}
   * @name Numerical settings for the curl-curl term
   * @{
   *
   * \var curlcurl_hodge
   * Set of parameters for the discrete Hodge operator related to the
   * curl-curl operator
   *
   * \var curlcurl_property
   * Pointer to the property related to the curl-curl term
   */

  cs_param_hodge_t              curlcurl_hodge;
  cs_property_t                *curlcurl_property;

  /*!
   * @}
   * @name Numerical settings for the grad-div term
   * @{
   *
   * \var graddiv_hodge
   * Set of parameters for the discrete Hodge operator related to the grad-div
   * operator
   *
   * \var graddiv_property
   * Pointer to the property related to the grad-div term
   */

  cs_param_hodge_t              graddiv_hodge;
  cs_property_t                *graddiv_property;

  /*!
   * @}
   * @name Numerical settings for the advection term
   * @{
   *
   * \var adv_formulation
   * Type of formulation (conservative, non-conservative...) for the advective
   * term
   *
   * \var adv_scheme
   * Numerical scheme used for the discretization of the advection term
   *
   * \var upwind_portion
   * Value between 0. and 1. (0: centered scheme, 1: pure upwind scheme)
   * Introduce a constant portion of upwinding in a centered scheme
   * Only useful if the advection scheme is set to
   *
   * \var adv_field
   * Pointer to the \ref cs_adv_field_t structure associated to the advection
   * term
   */

  cs_param_advection_form_t     adv_formulation;
  cs_param_advection_scheme_t   adv_scheme;
  cs_real_t                     upwind_portion;
  cs_adv_field_t               *adv_field;

  /*!
   * @}
   * @name Numerical settings for the reaction term
   * The contribution of a reaction term to the algebraic system is either at
   * the left-hand and/or right-hand side according to the choice of time
   * scheme
   * @{
   *
   * \var reaction_hodge
   * Set of parameters for the discrete Hodge operator related to the reaction
   * terms
   *
   * \var n_reaction_terms
   * Number of reaction terms to consider.
   *
   * \var reaction_properties
   * List of properties associated to each reaction term
   */

  cs_param_hodge_t              reaction_hodge;
  int                           n_reaction_terms;
  cs_property_t               **reaction_properties;

  /*!
   * @}
   * @name Definition of the related source terms
   * The contribution of a source term to the algebraic system is always at
   * right-hand side whatever is the choice of time scheme
   * @{
   *
   * \var n_source_terms
   * Number of source terms to consider.
   *
   * \var source_terms
   * List of definition of each source term
   */

  int                           n_source_terms;
  cs_xdef_t                   **source_terms;

  /*!
   * @}
   * @name Enforcement of values inside the computational domain
   *
   * This is different from the enforcement of boundary conditions but rely on
   * the same algebraic manipulation.
   * Two mechanisms are possible.
   *
   * 1) CELL SELECTION: defined a selection of cells and then
   * automatically built the related selection of degrees of freedom and
   * assigned a value to each selected degrees of freedom
   *
   * 2) DOF SELECTION: defined a selection of degrees of freedom (DoFs) and
   *    assign a values to a selection of degrees of freedom inside the domain
   *
   * @{
   *
   * \var enforcement_type
   * Flag specifying which kind of enforcement to perform
   *
   * \var enforcement_ref_value
   * Reference value to use. Avod to allocate an array with the same value
   * for all selected entities
   *
   * \var n_enforced_cells
   * Number of selected cells related to an enforcement
   *
   * \var enforced_cell_ids
   * List of selected cell ids related to an enforcement
   *
   * \var n_enforced_dofs
   * Number of degrees of freedom (DoFs) to enforce
   *
   * \var enforced_dof_ids
   * List of related DoF ids
   *
   * \var enforced_dof_values
   * List of related values to enforce
   */

  cs_flag_t                   enforcement_type;
  cs_real_t                  *enforcement_ref_value;

  cs_lnum_t                   n_enforced_cells;
  cs_lnum_t                  *enforced_cell_ids;
  cs_real_t                  *enforced_cell_values;

  cs_lnum_t                   n_enforced_dofs;
  cs_lnum_t                  *enforced_dof_ids;
  cs_real_t                  *enforced_dof_values;

  /*!
   * @}
   * @name Settings related to the resolution of the algebraic system
   * @{
   *
   * \var sles_param
   * Set of parameters to specify how to to solve the algebraic
   * - iterative solver
   * - preconditioner
   * - tolerance...
   */

  cs_param_sles_t              sles_param;

  /*!
   * @}
   * @name Settings related to the optimization of the performance
   * @{
   *
   * \var omp_assembly_choice
   * When OpenMP is active, choice of parallel reduction for the assembly
   */

  cs_param_assemble_omp_strategy_t     omp_assembly_choice;

  /*! @} */

} cs_equation_param_t;

/*! \enum cs_equation_key_t
 *  \brief List of available keys for setting the parameters of an equation
 *
 * \var CS_EQKEY_ADV_FORMULATION
 * Kind of formulation of the advective term. Available choices are:
 * - "conservative"
 * - "non_conservative"
 *
 * \var CS_EQKEY_ADV_SCHEME
 * Type of numerical scheme for the advective term. The available choices
 * depend on the space discretization scheme.
 * - "upwind" (cf. \ref CS_PARAM_ADVECTION_SCHEME_UPWIND)
 * - "centered" (cf. \ref CS_PARAM_ADVECTION_SCHEME_CENTERED)
 * - "mix_centered_upwind" (\ref CS_PARAM_ADVECTION_SCHEME_MIX_CENTERED_UPWIND)
 * - "samarskii" --> switch smoothly between an upwind and a centered scheme
 *   thanks to a weight depending on the Peclet number. (cf.
 * \ref CS_PARAM_ADVECTION_SCHEME_SAMARSKII). Only for CDO-Vb schemes.
 * - "sg" --> closely related to "samarskii" but with a different definition of
 *   the weight (cf. \ref CS_PARAM_ADVECTION_SCHEME_SG). Only for CDO-Vb schemes
 * - "cip" --> means "continuous interior penalty" (only for CDOVCB schemes).
 *   Enable a better accuracy. (cf. \ref CS_PARAM_ADVECTION_SCHEME_CIP)
 * - "cip_cw" --> means "continuous interior penalty" (only for CDOVCB schemes).
 *   Enable a better accuracy.
 *   Consider a cellwise approximation of the advection field.
 *   (cf. \ref CS_PARAM_ADVECTION_SCHEME_CIP_CW)
 *
 * \var CS_EQKEY_ADV_UPWIND_PORTION
 * Value between 0 and 1 specifying the portion of upwind added to a centered
 * discretization.
 *
 * \var CS_EQKEY_AMG_TYPE
 * Specify which type of algebraic multigrid (AMG) to choose.
 * Available choices are:
 * - "none" --> (default) No predefined AMG solver
 * - "boomer" --> Boomer AMG multigrid from the Hypre library
 * - "gamg" --> GAMG multigrid from the PETSc library
 * - "pcmg" --> MG multigrid as preconditioner from the PETSc library
 * - "v_cycle" --> Code_Saturne's in house multigrid with a V-cycle strategy
 * - "k_cycle" --> Code_Saturne's in house multigrid with a K-cycle strategy
 * WARNING: For "boomer" and "gamg",one needs to install Code_Saturne with
 * PETSc in this case
 *
 * \var CS_EQKEY_BC_ENFORCEMENT
 * Set the type of enforcement of the boundary conditions.
 * Available choices are:
 * - "algebraic": Modify the linear system so as to add the contribution of the
 * Dirichlet in the right-hand side and replace the line associated to a
 * Dirichlet BC by identity. This is a good choice for pure diffusinon or pure
 * convection problem.
 * - "penalization": Add a huge penalization coefficient on the diagonal term
 * of the line related to DoFs associated a Dirichlet BC. The right-hand side is
 * also modified to take into account this penalization. Be aware that it may
 * worsen the matrix conditioning.
 * - "weak": weak enforcement using the Nitsche method (there is also
 * penalization term but the scaling is such that a moderate value (1-100) of
 * the penalization coefficient is sufficient). This a good choice for
 * convection/diffusion problem.
 * - "weak_sym": Same as the "weak" option but with a modification so as to
 * add a symmetric contribution to the system. If the problem yields a symmetric
 * matrix. This choice is more relevant than "weak". This a good choice for a
 * diffusion problem.
 *
 * For HHO and CDO-Face based schemes, only the "penalization" and "algebraic"
 * technique is available up to now.
 *
 * \var CS_EQKEY_BC_QUADRATURE
 * Set the quadrature algorithm used for evaluating integral quantities on
 * faces or volumes. Available choices are:
 * - "bary"    used the barycenter approximation
 * - "higher"  used 4 Gauss points for approximating the integral
 * - "highest" used 5 Gauss points for approximating the integral
 *
 * Remark: "higher" and "highest" implies automatically a subdivision into
 * tetrahedra of each cell.
 *
 * \var CS_EQKEY_BC_STRONG_PENA_COEFF
 * Set the value of the penalization coefficient when "penalization" is
 * activated The default value is 1e12.
 * cf. \ref CS_PARAM_BC_ENFORCE_PENALIZED
 *
 * \var CS_EQKEY_BC_WEAK_PENA_COEFF
 * Set the value of the penalization coefficient when "weak" or "weak_sym" is
 * activated. The default value is 100.
 * cf. \ref CS_PARAM_BC_ENFORCE_WEAK_NITSCHE
 * or  \ref CS_PARAM_BC_ENFORCE_WEAK_SYM
 *
 * \var CS_EQKEY_DOF_REDUCTION
 * Set how is defined each degree of freedom (DoF).
 * - "de_rham" (default): Evaluation at vertices for potentials, integral
 *   along a line for circulations, integral across the normal component of a
 *   face for fluxes and integral inside a cell for densities
 * - "average": DoF are defined as the average on the element
 *
 * \var CS_EQKEY_EXTRA_OP
 * Set the additional post-processing to perform. Available choices are:
 * - "balance"  post-process the balance result in each control volume for
 *              each main term of an equation (diffusion, convection, time...)
 * - "peclet"  post-process an estimation of the Peclet number in each cell
 * - "upwind_coef"  post-process an estimation of the upwinding coefficient
 * - "normal_flux"  post-process the normal flux across boundary faces
 *
 * \var CS_EQKEY_HODGE_DIFF_ALGO
 * Set the algorithm used for building the discrete Hodge operator used
 * in the diffusion term. Available choices are:
 * - "voronoi" --> leads to a diagonal discrete Hodge operator but is not
 *   consistent for all meshes. Require an "orthogonal" (or admissible) mesh;
 * - "cost" --> (default for diffusion) is more robust (i.e. it handles more
 *   general meshes but is is less efficient)
 * - "wbs" --> is robust and accurate but is limited to the reconstruction of
 *   potential-like degrees of freedom and needs a correct computation of the
 *   cell barycenter
 *
 * \var CS_EQKEY_HODGE_DIFF_COEF
 * This key is only useful if CS_EQKEY_HODGE_{TIME, DIFF, REAC}_ALGO is set to
 * "cost".
 * keyval is either a name or a value:
 * - "dga" corresponds to the value \f$ 1./3. \f$
 * - "sushi" corresponds to the value \f$1./\sqrt(3.)\f$
 * - "gcr"  corresponds to the value \f$1\f$.
 * - or "1.5", "9" for instance
 *
 * \var CS_EQKEY_HODGE_TIME_ALGO
 * Set the algorithm used for building the discrete Hodge operator used
 * in the unsteady term. Available choices are:
 * - "voronoi" --> (default) leads to diagonal discrete Hodge operator. It acts
 *   as a mass lumping.
 * - "wbs" is robust and accurate but is less efficient. It needs a correct
 *   computation of the cell barycenter
 *
 * \var CS_EQKEY_HODGE_REAC_ALGO
 * Set the algorithm used for building the discrete Hodge operator used in the
 * reaction term. Available choices are:
 * - "voronoi" --> leads to diagonal discrete Hodge operator (similar to a
 * lumping).
 * - "wbs" --> (default) is robust and accurate but is limited to the
 * reconstruction of potential-like degrees of freedom and needs a correct
 * computation of the cell barycenter
 *
 * \var CS_EQKEY_ITSOL
 * Specify the iterative solver for solving the linear system related to an
 * equation. Avalaible choices are:
 * - "jacobi"           --> simpliest iterative solver
 * - "gauss_seidel"     --> Gauss-Seidel algorithm
 * - "sym_gauss_seidel" --> Symmetric version of Gauss-Seidel algorithm;
 *                          one backward and forward sweep
 * - "cg"               --> conjuguate gradient algorithm
 * - "fcg"              --> flexible version of the conjuguate gradient
 *                          algorithm used when the preconditioner can change
 *                          iteration by iteration
 * - "bicg"             --> Bi-CG algorithm (for non-symmetric linear systems)
 * - "bicgstab2"        --> BiCG-Stab2 algorithm (for non-symmetric linear
 *                          systems)
 * - "cr3"              --> a 3-layer conjugate residual solver (when "cs" is
 *                          chosen as the solver family)
 * - "gmres"            --> robust iterative solver. Not the best choice if the
 *                          system is easy to solve
 * - "amg"              --> algebraic multigrid iterative solver. Good choice
 *                          for a scalable solver related to symmetric positive
 *                          definite system.
 * - "minres"           --> Solver of choice for symmetric indefinite systems
 * - "mumps"            --> Direct solver (very robust but memory consumming)
 *                          via PETSc only. LU factorization.
 * - "mumps_ldlt"       --> Direct solver (very robust but memory consumming)
 *                          via PETSc only. LDLT factorization.
 *
 * \var CS_EQKEY_ITSOL_EPS
 * Tolerance factor for stopping the iterative processus for solving the
 * linear system related to an equation
 * - Example: "1e-10"
 *
 * \var CS_EQKEY_ITSOL_MAX_ITER
 * Maximum number of iterations for solving the linear system
 * - Example: "2000"
 *
 * \var CS_EQKEY_ITSOL_RESNORM_TYPE
 * Normalized or not the residual before testing if one continues iterating
 * for solving the linear system. This normalization is performed before
 * applying the boundary conditions to avoid handling the penalization of
 * Dirichlet boundary conditions. If the RHS norm is equal to zero, then
 * the "vol_tot" option is used for rescaling the residual.
 *
 * Available choices are:
 * "false" or "none" (default)
 * "rhs"
 * "weighted_rhs" or "weighted"
 * "filtered_rhs" or "fieltered_rhs"
 *
 * \var CS_EQKEY_OMP_ASSEMBLY_STRATEGY
 * Choice of the way to perform the assembly when OpenMP is active
 * Available choices are:
 * - "atomic" or "critical"
 *
 * \var CS_EQKEY_PRECOND
 * Specify the preconditioner associated to an iterative solver. Available
 * choices are:
 * - "jacobi" --> diagonal preconditoner
 * - "block_jacobi" --> Only with PETSc
 * - "poly1" --> Neumann polynomial of order 1 (only with Code_Saturne)
 * - "poly2" --> Neumann polynomial of order 2 (only with Code_Saturne)
 * - "ssor" --> symmetric successive over-relaxation (only with PETSC)
 * - "ilu0" --> incomplete LU factorization (only with PETSc)
 * - "icc0" --> incomplete Cholesky factorization (for symmetric matrices and
 *   only with PETSc)
 * - "amg" --> algebraic multigrid
 * - "amg_block" --> algebraic multigrid by block (useful for vector-valued
 *                   equations)
 *
 * \var CS_EQKEY_SLES_VERBOSITY
 * Level of details written by the code for the resolution of the linear system
 * - Examples: "0", "1", "2" or higher
 *
 * \var CS_EQKEY_SPACE_SCHEME
 * Set the space discretization scheme. Available choices are:
 * - "cdo_vb"  for CDO vertex-based scheme
 * - "cdo_vcb" for CDO vertex+cell-based scheme
 * - "cdo_fb"  for CDO face-based scheme
 * - "cdo_eb"  for CDO edge-based scheme
 * - "hho_p1"  for HHO schemes with \f$\mathbb{P}_1\f$ polynomial approximation
 * - "hho_p2"  for HHO schemes with \f$\mathbb{P}_2\f$ polynomial approximation
 *
 * \var CS_EQKEY_TIME_SCHEME
 * Set the scheme for the temporal discretization. Available choices are:
 * - "euler_implicit": first-order in time (inconditionnally stable)
 * - "euler_explicit":
 * - "crank_nicolson": second_order in time
 * - "theta_scheme": generic time scheme. One recovers "euler_implicit" with
 *   theta equal to "1", "explicit" with "0", "crank_nicolson" with "0.5"
 *
 * \var CS_EQKEY_TIME_THETA
 * Set the value of theta. Only useful if CS_EQKEY_TIME_SCHEME is set to
 * "theta_scheme"
 * - Example: "0.75" (keyval must be between 0 and 1)
 *
 * \var CS_EQKEY_VERBOSITY
 * Set the level of details written by the code for an equation.
 * The higher the more detailed information is displayed.
 * - "0" (default)
 * - "1" detailed setup resume and coarse grain timer stats
 * - "2" fine grain for timer stats
 *
 */

typedef enum {

  CS_EQKEY_ADV_FORMULATION,
  CS_EQKEY_ADV_SCHEME,
  CS_EQKEY_ADV_UPWIND_PORTION,
  CS_EQKEY_AMG_TYPE,
  CS_EQKEY_BC_ENFORCEMENT,
  CS_EQKEY_BC_QUADRATURE,
  CS_EQKEY_BC_STRONG_PENA_COEFF,
  CS_EQKEY_BC_WEAK_PENA_COEFF,
  CS_EQKEY_DO_LUMPING,
  CS_EQKEY_DOF_REDUCTION,
  CS_EQKEY_EXTRA_OP,
  CS_EQKEY_HODGE_DIFF_ALGO,
  CS_EQKEY_HODGE_DIFF_COEF,
  CS_EQKEY_HODGE_TIME_ALGO,
  CS_EQKEY_HODGE_REAC_ALGO,
  CS_EQKEY_ITSOL,
  CS_EQKEY_ITSOL_EPS,
  CS_EQKEY_ITSOL_MAX_ITER,
  CS_EQKEY_ITSOL_RESNORM_TYPE,
  CS_EQKEY_OMP_ASSEMBLY_STRATEGY,
  CS_EQKEY_PRECOND,
  CS_EQKEY_SLES_VERBOSITY,
  CS_EQKEY_SOLVER_FAMILY,
  CS_EQKEY_SPACE_SCHEME,
  CS_EQKEY_TIME_SCHEME,
  CS_EQKEY_TIME_THETA,
  CS_EQKEY_VERBOSITY,

  CS_EQKEY_N_KEYS

} cs_equation_key_t;

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the flag related to a cs_equation_param_t structure
 *
 * \param[in, out] eqp          pointer to a \ref cs_equation_param_t
 * \param[in]      flag         flag to add
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_equation_param_set_flag(cs_equation_param_t     *eqp,
                           cs_flag_t                flag)
{
  assert(eqp != NULL);
  eqp->flag |= flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs a diffusion term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_diffusion(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_DIFFUSION)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs a curl-curl term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_curlcurl(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_CURLCURL)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs a grad-div term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_graddiv(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_GRADDIV)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs a convection term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_convection(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_CONVECTION)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs a reaction term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_reaction(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_REACTION)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs an unsteady term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_time(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation needs a source term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_sourceterm(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->n_source_terms > 0)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation has an internal enforcement
 *         of the degrees of freedom
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_internal_enforcement(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_FORCE_VALUES)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a \ref cs_equation_param_t structure has its name member
 *         equal to the given name
 *
 * \param[in] eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in] name       name of the equation
 *
 * \return true if the given eqp has the same name the one given as parameter
 *         otherwise false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_name(cs_equation_param_t   *eqp,
                           const char            *name)
{
  if (eqp == NULL)
    return false;
  if (eqp->name == NULL)
    return false;
  if (strcmp(eqp->name, name) == 0)
    return true;
  else
    return false;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a \ref cs_equation_param_t structure
 *
 * \param[in] name          name of the equation
 * \param[in] type          type of equation
 * \param[in] dim           dim of the variable associated to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return a pointer to a new allocated \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_create_param(const char            *name,
                         cs_equation_type_t     type,
                         int                    dim,
                         cs_param_bc_type_t     default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the settings from one \ref cs_equation_param_t structure to
 *         another one
 *
 * \param[in]      ref   pointer to the reference \ref cs_equation_param_t
 * \param[in, out] dst   pointer to the \ref cs_equation_param_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_update_from(const cs_equation_param_t   *ref,
                              cs_equation_param_t         *dst);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_equation_param_t
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_free_param(cs_equation_param_t     *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a \ref cs_equation_param_t
 *         structure
 *
 * \param[in, out]  eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_param(cs_equation_param_t   *eqp,
                      cs_equation_key_t      key,
                      const char            *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Settings are related to this equation.
 *
 * \param[in, out]  eqp       pointer to a \ref cs_equation_param_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_sles(cs_equation_param_t     *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last modification of the cs_equation_param_t structure before
 *         launching the computation
 *
 * \param[in, out]  eqp      pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_last_stage(cs_equation_param_t   *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a \ref cs_equation_param_t structure
 *
 * \param[in]  eqp      pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_summary_param(const cs_equation_param_t  *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here a constant value is set to all the entities belonging to the
 *         given mesh location
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]       val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_value(cs_equation_param_t    *eqp,
                            const char             *z_name,
                            cs_real_t              *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the value related to all the entities belonging to the
 *         given mesh location is such that the integral over these cells
 *         returns the requested quantity
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]       quantity  quantity to distribute over the mesh location
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_qov(cs_equation_param_t    *eqp,
                          const char             *z_name,
                          double                  quantity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this
 *         equation. This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_analytic(cs_equation_param_t    *eqp,
                               const char             *z_name,
                               cs_analytic_func_t     *analytic,
                               void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a boundary condition from an existing \ref cs_xdef_t structure
 *         The lifecycle of the cs_xdef_t structure is now managed by the
 *         current \ref cs_equation_param_t structure.
 *
 * \param[in, out] eqp    pointer to a cs_equation_param_t structure
 * \param[in]      xdef   pointer to the \ref cs_xdef_t structure to transfer
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_add_xdef_bc(cs_equation_param_t        *eqp,
                        cs_xdef_t                  *xdef);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       values    pointer to a array storing the values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_value(cs_equation_param_t         *eqp,
                            const cs_param_bc_type_t     bc_type,
                            const char                  *z_name,
                            cs_real_t                   *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                            (true or false)
 * \param[in]       index     optional pointer to the array index
 *
 * \return a pointer to the new allocated \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_array(cs_equation_param_t        *eqp,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_flag_t                   loc,
                            cs_real_t                  *array,
                            bool                        is_owner,
                            cs_lnum_t                  *index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation param structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function defining the value
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_analytic(cs_equation_param_t        *eqp,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *z_name,
                               cs_analytic_func_t         *analytic,
                               void                       *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a sliding boundary
 *         condition related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the related boundary zone
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_sliding_condition(cs_equation_param_t     *eqp,
                                  const char              *z_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the Laplacian operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *         Laplacian means div-grad (either for vector-valued or scalar-valued
 *         equations)
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_diffusion(cs_equation_param_t   *eqp,
                          cs_property_t         *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the curl-curl operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_curlcurl(cs_equation_param_t   *eqp,
                         cs_property_t         *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the grad-div operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_graddiv(cs_equation_param_t   *eqp,
                        cs_property_t         *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the time derivative operator for
 *         the equation associated to the given \ref cs_equation_param_t
 *         structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_time(cs_equation_param_t   *eqp,
                     cs_property_t         *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the advection operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      adv_field  pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_advection(cs_equation_param_t   *eqp,
                          cs_adv_field_t        *adv_field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the reaction operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 *
 * \return the id related to the reaction term
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_add_reaction(cs_equation_param_t   *eqp,
                         cs_property_t         *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by a constant value
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or
 *                           "" all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_val(cs_equation_param_t    *eqp,
                                   const char             *z_name,
                                   cs_real_t              *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by an analytical function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      func      pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_analytic(cs_equation_param_t    *eqp,
                                        const char             *z_name,
                                        cs_analytic_func_t     *func,
                                        void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by a DoF function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  location of element ids given as parameter
 * \param[in]      func      pointer to a DoF function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_dof_func(cs_equation_param_t    *eqp,
                                        const char             *z_name,
                                        cs_flag_t               loc_flag,
                                        cs_dof_func_t          *func,
                                        void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by an array.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                           (true or false)
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_array(cs_equation_param_t    *eqp,
                                     const char             *z_name,
                                     cs_flag_t               loc,
                                     cs_real_t              *array,
                                     bool                    is_owner,
                                     cs_lnum_t              *index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value of degrees of freedom located at
 *         mesh vertices.
 *         The spatial discretization scheme for the given equation has to be
 *         CDO-Vertex based or CDO-Vertex+Cell-based schemes.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp         pointer to a cs_equation_param_t structure
 * \param[in]      n_elts      number of vertices to enforce
 * \param[in]      elt_ids     list of vertices
 * \param[in]      ref_value   ignored if NULL
 * \param[in]      elt_values  list of associated values, ignored if NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforce_vertex_dofs(cs_equation_param_t    *eqp,
                                cs_lnum_t               n_elts,
                                const cs_lnum_t         elt_ids[],
                                const cs_real_t         ref_value[],
                                const cs_real_t         elt_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value related to the degrees of freedom
 *         associated to the list of selected cells.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp         pointer to a cs_equation_param_t structure
 * \param[in]      n_elts      number of selected cells
 * \param[in]      elt_ids     list of cell ids
 * \param[in]      ref_value   ignored if NULL
 * \param[in]      elt_values  list of associated values, ignored if NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforce_by_cell_selection(cs_equation_param_t    *eqp,
                                      cs_lnum_t               n_elts,
                                      const cs_lnum_t         elt_ids[],
                                      const cs_real_t         ref_value[],
                                      const cs_real_t         elt_values[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_PARAM_H__ */
