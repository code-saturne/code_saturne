#ifndef __CS_EQUATION_PARAM_H__
#define __CS_EQUATION_PARAM_H__

/*============================================================================
 * Header to handle specific settings related to a cs_equation_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_property.h"
#include "cs_advection_field.h"
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
 * \brief Diffusion term is needed
 *
 * \def CS_EQUATION_REACTION
 * \brief Reaction term is needed
 *
 */
#define CS_EQUATION_LOCKED        (1 <<  0)  //  1
#define CS_EQUATION_UNSTEADY      (1 <<  1)  //  2
#define CS_EQUATION_CONVECTION    (1 <<  2)  //  4
#define CS_EQUATION_DIFFUSION     (1 <<  3)  //  8
#define CS_EQUATION_REACTION      (1 <<  4)  // 16

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
 */

#define CS_EQUATION_POST_BALANCE     (1 << 0) /* 1 */
#define CS_EQUATION_POST_PECLET      (1 << 1) /* 2 */
#define CS_EQUATION_POST_UPWIND_COEF (1 << 2) /* 4 */

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
 * \var CS_EQUATION_TYPE_GROUNDWATER,
 * Equation related to the groundwater flow module
 *
 * \var CS_EQUATION_TYPE_PREDEFINED
 * Predefined equation (most part of the setting is already done)
 * - Example: equation for the wall distance
 */

typedef enum {

  CS_EQUATION_TYPE_USER,
  CS_EQUATION_TYPE_GROUNDWATER,
  CS_EQUATION_TYPE_PREDEFINED,
  CS_EQUATION_N_TYPES

} cs_equation_type_t;

/*! \enum cs_equation_solver_class_t
 *  \brief Class of iterative solvers to consider for solver the linear system
 *
 * \var CS_EQUATION_SOLVER_CLASS_CS
 * Iterative solvers available in Code_Saturne
 *
 * \var CS_EQUATION_SOLVER_CLASS_PETSC
 * Solvers available in Code_Saturne
 *
 * \var CS_EQUATION_N_SOLVER_CLASSES
 */

typedef enum {

  CS_EQUATION_SOLVER_CLASS_CS,
  CS_EQUATION_SOLVER_CLASS_PETSC,
  CS_EQUATION_N_SOLVER_CLASSES

} cs_equation_solver_class_t;

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
  int                  sles_verbosity; /*!< Level of detail for SLES output */

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
   * \var enforcement
   * Type of enforcement for the Dirichlet boundary conditions.
   * See \ref cs_param_bc_enforce_t for more details.
   *
   * \var n_bc_defs
   * Number of boundary conditions which are defined for this equation
   *
   * \var bc_defs
   * Pointers to the definitions of the boundary conditions
   */

  cs_param_bc_type_t            default_bc;
  cs_param_bc_enforce_t         enforcement;
  int                           n_bc_defs;
  cs_xdef_t                   **bc_defs;

  /*!
   * @}
   * @name Settings related to the resolution of the algebraic system
   * @{
   *
   * \var solver_class
   * Class of solver available to solve the algebraic system
   *
   * \var itsol_info
   * Set of parameters to specify how to to solve the algebraic
   * - iterative solver
   * - preconditionner
   * - tolerance...
   */

  cs_equation_solver_class_t    solver_class;
  cs_param_itsol_t              itsol_info;

  /*!
   * @}
   * @name Numerical settings for the time-dependent parameters
   * @{
   *
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
   * \var do_lumping
   * Perform a mass lumping on the matrix related to the time discretization.
   * - true or false
   */

  cs_param_hodge_t              time_hodge;
  cs_property_t                *time_property;
  cs_param_time_scheme_t        time_scheme;
  cs_real_t                     theta;
  bool                          do_lumping;

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
   * @name Numerical settings for the diffusion term
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
   * \var adv_field
   * Pointer to the \cs_adv_field_t structure associated to the advection term
   */

  cs_param_advection_form_t     adv_formulation;
  cs_param_advection_scheme_t   adv_scheme;
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

  /*! @} */

} cs_equation_param_t;

/*! \enum cs_equation_key_t
 *  \brief List of available keys for setting the parameters of an equation
 *
 * \var CS_EQKEY_SPACE_SCHEME
 * Set the space discretization scheme. Available choices are:
 * - "cdo_vb"  for CDO vertex-based scheme
 * - "cdo_vcb" for CDO vertex+cell-based scheme
 * - "cdo_fb"  for CDO face-based scheme
 * - "hho_p1"  for HHO schemes with \f$\mathbb{P}_1\f$ polynomial approximation
 * - "hho_p2"  for HHO schemes with \f$\mathbb{P}_2\f$ polynomial approximation
 *
 * \var CS_EQKEY_DOF_REDUCTION
 * Set how is defined each degree of freedom (DoF).
 * - "de_rham" (default): Evaluation at vertices for potentials, integral
 *   along a line for circulations, integral across the normal component of a
 *   face for fluxes and integral inside a cell for densities
 * - "average": DoF are defined as the average on the element
 *
 * \var CS_EQKEY_VERBOSITY
 * Set the level of details written by the code for an equation.
 * The higher the more detailed information is displayed.
 * - "0" (default)
 * - "1" detailed setup resume and coarse grain timer stats
 * - "2" fine grain for timer stats
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
 * \var CS_EQKEY_HODGE_TIME_ALGO
 * Set the algorithm used for building the discrete Hodge operator used
 * in the unsteady term. Available choices are:
 * - "voronoi" --> (default) leads to diagonal discrete Hodge operator. It acts
 *   as a mass lumping.
 * - "wbs" is robust and accurate but is less efficient. It needs a correct
 *   computation of the cell barycenter
 *
 * \var CS_EQKEY_HODGE_REAC_ALGO
 * Set the algorithm used for building the discrete Hodge operator used
 * in the reaction term. Available choices are:
 * - "voronoi" --> leads to diagonal discrete Hodge operator but is not
 *   consistent for all meshes.  Require an "orthogonal" (or admissible) mesh;
 * - "wbs" --> (default) is robust and accurate but is limited to the
 *   reconstruction of potential-like degrees of freedom and needs a correct
 *   computation of the cell barycenter
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
 * \var CS_EQKEY_HODGE_TIME_COEF
 * This key is only useful if CS_EQKEY_HODGE_{TIME, DIFF, REAC}_ALGO is set to
 * "cost".
 * keyval is either a name or a value:
 * - "dga" corresponds to the value \f$ 1./3. \f$
 * - "sushi" corresponds to the value \f$1./\sqrt(3.)\f$
 * - "gcr"  corresponds to the value \f$1\f$.
 * - or "1.5", "9" for instance
 *
 * \var CS_EQKEY_HODGE_REAC_COEF
 * This key is only useful if CS_EQKEY_HODGE_{TIME, DIFF, REAC}_ALGO is set to
 * "cost".
 * keyval is either a name or a value:
 * - "dga" corresponds to the value \f$ 1./3. \f$
 * - "sushi" corresponds to the value \f$1./\sqrt(3.)\f$
 * - "gcr"  corresponds to the value \f$1\f$.
 * - or "1.5", "9" for instance
 *
 * \var CS_EQKEY_SOLVER_FAMILY
 * Specify which class of solver are possible. Available choices are:
 * - "none" --> (default) List of possible iterative solvers are those of
 *   Code_Saturne,
 * - "petsc" --> List of possible iterative solvers are those of the PETSc
 *   library. WARNING: one needs to install Code_Saturne with PETSc in this case
 *
 * \var CS_EQKEY_AMG_TYPE
 * Specify which type of algebraic multigrid (AMG) to choose.
 * Available choices are:
 * - "none" --> (default) No predefined AMG solver
 * - "boomer" --> Boomer AMG multigrid from the Hypre library
 * - "gamg" --> GAMG multigrid from the PETSc library
 * - "v_cycle" --> Code_Saturne's in house multigrid with a V-cycle strategy
 * - "k_cycle" --> Code_Saturne's in house multigrid with a K-cycle strategy
 * WARNING: For "boomer" and "gamg",one needs to install Code_Saturne with
 * PETSc in this case
 *
 * \var CS_EQKEY_ITSOL
 * Specify the iterative solver for solving the linear system related to an
 * equation. Avalaible choices are:
 * - "jacobi" --> simpliest algorithm
 * - "gauss_seidel" --> Gauss-Seidel algorithm
 * - "sym_gauss_seidel" --> Symmetric version of Gauss-Seidel algorithm;
 *                          one backward and forward sweep
 * - "cg" --> (default) the standard conjuguate gradient algorithm
 * - "fcg" --> flexible version of the conjuguate gradient algorithm used
 *             when the preconditioner can change iteration by iteration
 * - "bicg" --> Bi-CG algorithm (for non-symmetric linear systems)
 * - "bicgstab2" --> BiCG-Stab2 algorithm (for non-symmetric linear systems)
 * - "cr3" --> a 3-layer conjugate residual solver (when "cs" is chosen as the
 *    solver family)
 * - "gmres" --> a robust iterative solver but slower as previous one if the
 *   system is not difficult to solve
 * - "amg" --> an algebraic multigrid iterative solver. Good choice for a
 * symmetric positive definite system.
 *
 * \var CS_EQKEY_PRECOND
 * Specify the preconditionner associated to an iterative solver. Available
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
 * \var CS_EQKEY_ITSOL_RESNORM
 * Normalized or not the residual before testing if one continues iterating
 * for solving the linear system. Available choices are:
 * - "true" or "false"
 *
 * \var CS_EQKEY_SLES_VERBOSITY
 * Level of details written by the code for the resolution of the linear system
 * - Examples: "0", "1", "2" or higher
 *
 * \var CS_EQKEY_BC_ENFORCEMENT
 * Set the type of enforcement of the boundary conditions.
 * Available choices are:
 * - "penalization" --> weak enforcement using a huge penalization coefficient
 * - "weak" --> weak enforcement using the Nitsche method
 * - "weak_sym" --> weak enforcement keeping the symmetry of the system
 *
 * For HHO and CDO-Face based schemes, only the penalization technique is
 * available.
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
 * \var CS_EQKEY_TIME_SCHEME
 * Set the scheme for the temporal discretization. Available choices are:
 * - "implicit": first-order in time (inconditionnally stable)
 * - "explicit":
 * - "crank_nicolson": second_order in time
 * - "theta_scheme": generic time scheme. One recovers "implicit" with theta
 *   equal to "1", "explicit" with "0", "crank_nicolson" with "0.5"
 *
 * \var CS_EQKEY_TIME_THETA
 * Set the value of theta. Only useful if CS_EQKEY_TIME_SCHEME is set to
 * "theta_scheme"
 * - Example: "0.75" (keyval must be between 0 and 1)
 *
 * \var CS_EQKEY_ADV_FORMULATION
 * Kind of formulation of the advective term. Available choices are:
 * - "conservative"
 * - "non_conservative"
 *
 * \var CS_EQKEY_ADV_SCHEME
 * Type of numerical scheme for the advective term. The available choices
 * depend on the space discretization scheme.
 * - "upwind"
 * - "centered"
 * - "samarskii" --> switch smoothly betwwen an upwind and a centered scheme
 *   thanks to a weight depending on the Peclet number.
 * - "sg" --> closely related to "samarskii" but with a different definition of
 *   the weight.
 * - "cip" --> means "continuous interior penalty" (only for CDOVCB schemes).
 *   Enable a better accuracy.
 *
 * "sg" and "samarskii" are only available with CDOVB schemes
 *
 *
 * \var CS_EQKEY_EXTRA_OP
 * Set the additional post-processing to perform. Available choices are:
 * - "peclet" --> post-process an estimation of the Peclet number in each cell
 * - "upwind_coef" --> post-process an estimation of the upwinding coefficient
 */

typedef enum {

  CS_EQKEY_ADV_FORMULATION,
  CS_EQKEY_ADV_SCHEME,
  CS_EQKEY_AMG_TYPE,
  CS_EQKEY_BC_ENFORCEMENT,
  CS_EQKEY_BC_QUADRATURE,
  CS_EQKEY_DOF_REDUCTION,
  CS_EQKEY_EXTRA_OP,
  CS_EQKEY_HODGE_DIFF_ALGO,
  CS_EQKEY_HODGE_DIFF_COEF,
  CS_EQKEY_HODGE_TIME_ALGO,
  CS_EQKEY_HODGE_TIME_COEF,
  CS_EQKEY_HODGE_REAC_ALGO,
  CS_EQKEY_HODGE_REAC_COEF,
  CS_EQKEY_ITSOL,
  CS_EQKEY_ITSOL_EPS,
  CS_EQKEY_ITSOL_MAX_ITER,
  CS_EQKEY_ITSOL_RESNORM,
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
  if (eqp == NULL)
    return false;
  if (eqp->flag & CS_EQUATION_DIFFUSION)
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
  if (eqp == NULL)
    return false;
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
  if (eqp == NULL)
    return false;
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
  if (eqp == NULL)
    return false;
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
  if (eqp == NULL)
    return false;
  if (eqp->n_source_terms > 0)
    return true;
  else
    return false;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a \ref cs_equation_param_t structure has its name member
 *         equal to the given name
 *
 * \apram[in] eqp        pointer to a \ref cs_equation_param_t structure
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

 * \param[in]   eqp          pointer to a \ref cs_equation_param_t struct.
 * \param[in]   field_id     id of the cs_field_t struct. for this equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_sles(cs_equation_param_t     *eqp,
                           int                      field_id);

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
 * \param[in]       index     optional pointer to the array index
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_array(cs_equation_param_t        *eqp,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_flag_t                   loc,
                            cs_real_t                  *array,
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
 * \brief  Define and initialize a new structure to store parameters related
 *         to a diffusion term
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
 * \brief  Define and initialize a new structure to store parameters related
 *         to an unsteady term
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
 * \brief  Define and initialize a new structure to store parameters related
 *         to an advection term
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
 * \brief  Define and initialize a new structure to store parameters related
 *         to a reaction term
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
 * \brief  Define a new source term structure and initialize it by value
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or
 *                           "" all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_val(cs_equation_param_t    *eqp,
                                   const char             *z_name,
                                   cs_real_t              *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure and initialize it by an analytical
 *         function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_analytic(cs_equation_param_t    *eqp,
                                        const char             *z_name,
                                        cs_analytic_func_t     *ana,
                                        void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term defined by an array
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_array(cs_equation_param_t    *eqp,
                                     const char             *z_name,
                                     cs_flag_t               loc,
                                     cs_real_t              *array,
                                     cs_lnum_t              *index);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_PARAM_H__ */
