#ifndef __CS_EQUATION_PARAM_H__
#define __CS_EQUATION_PARAM_H__

/*============================================================================
 * Functions related to the structure cs_equation_param_t storing the settings
 * related to an equation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include "cs_enforcement.h"
#include "cs_hodge.h"
#include "cs_iter_algo.h"
#include "cs_param_cdo.h"
#include "cs_param_sles.h"
#include "cs_property.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @name Flags specifying the behavior of an equation (activated terms like the
 *       convection, diffusion, etc.).
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
 *        with div .grad
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
 *        interior degrees of freedom
 *
 * \def CS_EQUATION_INSIDE_SYSTEM
 * \brief The current equation settings belong to a system of equation
 *
 * \def CS_EQUATION_BUILD_HOOK
 * \brief Activate a build hook function to get a fine control of the
 *        discretization process during the cellwise building of the linear
 *        system.
 *        Need a function matching the cs_equation_build_hook_t prototype
 *
 * \def CS_EQUATION_USER_TRIGGERED
 * \brief The stage at which the equation is solved is driven by the user.  So,
 *        this user-defined equation is not solved at the same time as the
 *        other user-defined equations.
 */

#define CS_EQUATION_LOCKED         (1 <<  0)  /*    1 */
#define CS_EQUATION_UNSTEADY       (1 <<  1)  /*    2 */
#define CS_EQUATION_CONVECTION     (1 <<  2)  /*    4 */
#define CS_EQUATION_DIFFUSION      (1 <<  3)  /*    8 */
#define CS_EQUATION_CURLCURL       (1 <<  4)  /*   16 */
#define CS_EQUATION_GRADDIV        (1 <<  5)  /*   32 */
#define CS_EQUATION_REACTION       (1 <<  6)  /*   64 */
#define CS_EQUATION_FORCE_VALUES   (1 <<  7)  /*  128 */
#define CS_EQUATION_INSIDE_SYSTEM  (1 <<  8)  /*  256 */
#define CS_EQUATION_BUILD_HOOK     (1 <<  9)  /*  512 */
#define CS_EQUATION_USER_TRIGGERED (1 << 10)  /* 1024 */


/*!
 * @}
 * @name Flags specifying which extra operation related to an equation is
 *       needed.
 * @{
 *
 * \def CS_EQUATION_POST_BALANCE
 * \brief Compute and postprocess the equation balance
 *
 * \def CS_EQUATION_POST_PECLET
 * \brief Compute and postprocess the Peclet number
 *
 * \def CS_EQUATION_POST_UPWIND_COEF
 * \brief Postprocess the value of the upwinding coefficient. Only useful if
 *        variable coefficient controls the amount of upwinding. For instant,
 *        in schemes such as Schaffeter-Gummel, the upwind portion is driven by
 *        the local Peclet number.
 *
 * \def CS_EQUATION_POST_NORMAL_FLUX
 * \brief Postprocess the value of the normal flux across the boundary faces
 *
 */

#define CS_EQUATION_POST_BALANCE     (1 << 0) /* 1 */
#define CS_EQUATION_POST_PECLET      (1 << 1) /* 2 */
#define CS_EQUATION_POST_UPWIND_COEF (1 << 2) /* 4 */
#define CS_EQUATION_POST_NORMAL_FLUX (1 << 3) /* 8 */

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_equation_type_t
 *  \brief Type of equations managed by the solver
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
 * \var CS_EQUATION_TYPE_PREDEFINED
 * Predefined equation (most part of the setting is already done)
 * - Example: equation for the wall distance or ALE
 *
 * \var CS_EQUATION_TYPE_THERMAL
 * Equation related to the heat transfer
 *
 * \var CS_EQUATION_TYPE_SOLIDIFICATION
 * Equation related to the solidification module
 *
 * \var CS_EQUATION_TYPE_USER
 * User-defined equation (no interaction with energy and momentum equations)
 * Resolution is performed after all predefined equations have been solved.
 */

typedef enum {

  CS_EQUATION_TYPE_GROUNDWATER,
  CS_EQUATION_TYPE_MAXWELL,
  CS_EQUATION_TYPE_NAVSTO,
  CS_EQUATION_TYPE_PREDEFINED,
  CS_EQUATION_TYPE_THERMAL,
  CS_EQUATION_TYPE_SOLIDIFICATION,
  CS_EQUATION_TYPE_USER,

  CS_EQUATION_N_TYPES

} cs_equation_type_t;

/*----------------------------------------------------------------------------*/

/*! \struct cs_equation_param_t
 *  \brief Set of parameters to handle an unsteady convection-diffusion-reaction
 *         equation with term sources
 */

typedef struct {

  /*!
   * @name General settings
   * @{
   */

  char *restrict              name;    /*!< name of the equation */
  cs_equation_type_t          type;    /*!< type of equation: predefined... */
  int                         dim;     /*!< Dimension of the unknown */

  /*! \var verbosity
   *  Verbosity for the resolution (0 or 1 for a reasonable log size,  2 or more
   * for troubleshooting).
   *
   * \var iwarni
   * \deprecated use verbosity instead (iwarni is an alias to verbosity)
   */

  union {
    int                       verbosity;
    int                       iwarni;
  };

  /*! \var flag
   *  Flag to know if unsteady or diffusion or convection or reaction
   *  or source terms are activated or not
   */

  cs_flag_t                   flag;

  /*! \var post_flag
   *  Flag to determine if predefined post-treatments such as Peclet,
   *  are requested
   */

  cs_flag_t                   post_flag;

  /* Numerical settings */

  cs_param_space_scheme_t     space_scheme;  /*!< Space discretization scheme */
  cs_param_dof_reduction_t    dof_reduction; /*!< How is defined DoF */

  /*! \var space_poly_degree
   * Maximum degree of the polynomial basis
   */

  int                         space_poly_degree;

  /*!
   * @}
   * @name Legacy Settings
   * @{
   *
   * \var iconv

   * Indicate if the convection is taken into account (1) or not (0). By default,
   * 0 for the pressure or f in v2f model, 1 for the other unknowns.
   *
   * \var istat
   * Indicate whether unsteady terms are present (1) or not (0) in the matrices.
   * By default, 0 for the pressure or f in v2f model, 1 for the other unknowns.
   *
   * \var idircl
   * Indicate whether the diagonal of the matrix should be slightly shifted if
   * there is no Dirichlet boundary condition and if \ref istat = 0
   * (0: false / 1: true). Indeed, in such a case, the matrix for the general
   * advection/diffusion equation is singular. A slight shift in the diagonal
   * will make it invertible again.
   * By default, \ref idircl is set to 1 for all the unknowns, except
   * \f$\overline{f}\f$ in the v2f model (whose equation already contain another
   * diagonal term).
   * \remark this code is defined automatically based on the
   * presence of Dirichlet BCs.
   *
   * \var ndircl
   * Number of Dirichlet BCs
   *
   * \var idiff
   * Indifcate if diffusion is taken into account (1) or not (0).
   *
   * \var idifft
   * When diffusion is taken into account (\ref idiff = 1), indicate if the
   * turbulent diffusion is taken into account (\ref idifft = 1) or not (0).
   *
   * \var idften
   * Type of diffusivity flag (sum of mask constants defining if diffusivity is
   * isotropic, anisotropic, ... Masks are defined in \ref scalar_params).
   *
   * \var iswdyn
   * Dynamic relaxation type:
   * - 0 no dynamic relaxation
   * - 1 dynamic relaxation depending on \f$ \delta \varia^k \f$
   * - 2 dynamic relaxation depending on \f$ \delta \varia^k \f$
   *    and \f$ \delta \varia^{k-1} \f$.
   *
   * \var ischcv
   * Indicate the type of second-order convective scheme
   * - 0: Second Order Linear Upwind
   * - 1: Centered
   * - 2: Second Order with upwind-gradient reconstruction (SOLU)
   * - 3: Blending between SOLU and Centered scheme
   * - 4: NVD/TVD Scheme
   * Then "limiter_choice" keyword must be set:
   * - 0: Gamma
   * - 1: SMART
   * - 2: CUBISTA
   * - 3: SUPERBEE
   * - 4: MUSCL
   * - 5: MINMOD
   * - 6: CLAM
   * - 7: STOIC
   * - 8: OSHER
   * - 9: WASEB
   * --- VOF scheme ---
   * - 10: M-HRIC
   * - 11: M-CICSAM
   * - 12: M-STACS
   *
   * \var ibdtso
   * Backward differential scheme in time order.
   *
   * \var isstpc
   * Indicate whether a slope test should be used to switch from a second-order
   * to an upwind convective scheme under certain conditions, to ensure
   * stability.
   * - 0: slope test activated
   * - 1: slope test deactivated
   * - 2: continuous limiter ensuring boundedness (beta limiter)
   * The use of the slope test stabilises the calculation but may reduce the
   * spatial convergence order.
   *
   * \var nswrgr
   * Iteration limit for the iterative gradient reconstruction (\ref imrgra = 0).
   * If \ref imrgra = 0 and \ref nswrgr <= 1, gradients are not reconstructed.
   *
   * \var nswrsm
   * Iteration limit for the reconstruction of the right-hand sides of the
   * equations with a first-order scheme in time (standard case), the default
   * values are 2 for pressure and 1 for the other variables. With a
   * second-order scheme in time (\ref optcal::ischtp "ischtp" = 2) or LES, the
   * default values are 5 for pressure and 10 for the other variables.
   *
   * \var imvisf
   * Face viscosity field interpolation
   * - 0: aritmetic (default)
   * - 1: harmonic
   *
   * \var imrgra
   * Type of gradient reconstruction
   * - 0: iterative reconstruction of the non-orthogonalities
   * - 1: least squares method based on the first neighbor cells
   *      (those which share a face with the treated cell)
   * - 2, 3: least squares method using the extended neighborhood
   * - 4: Green-Gauss based using the least squares method
   *      (first neighbors) to compute face values
   * - 5, 6: Green-Gauss based using the least squares method with
   *      an extended neighborhood to compute face values\n
   * If the computation fails due to mesh quality aspects,
   * it is usually effective to use imrgra = 3, 5, or 6.
   *
   * \var imligr
   * Type of gradient limiter
   * - -1 (CS_GRADIENT_LIMIT_NONE): no limitation
   * - 0 (CS_GRADIENT_LIMIT_CELL): based on the neighbors
   * - 1 (CS_GRADIENT_LIMIT_FACE): superior order
   *
   * \var ircflu
   * Indicate whether the convective and diffusive fluxes at the faces should be
   * reconstructed:
   * - 0: no reconstruction
   * - 1: reconstruction
   * Deactivating the reconstruction of the fluxes can have a stabilizing
   * effect on the calculation. It is sometimes useful with the \f$ k-\epsilon
   * \f$ model, if the mesh is strongly non-orthogonal in the near-wall region,
   * where the gradients of k and \f$ \epsilon \f$ are strong. In such a case,
   * setting \ref ircflu = 0 will probably help (switching to a first order
   * convective scheme, \ref blencv = 0, for k and \f$ \epsilon \f$ might also
   * help in that case).
   *
   * \var iwgrec
   * Gradient calculation weighting
   * - 0: standard
   * - 1: weighted
   *
   * \var icoupl
   * Internal coupling indicator
   * - -1: not coupled (default)
   * -  1: coupled
   *
   * \var thetav
   * Value of \f$ \theta \f$ used to express at the second order the terms of
   * convection, diffusion and the source terms which are linear functions of
   * the solved variable, according to the formula
   * \f$ \phi^{n+\theta} = (1-\theta) \phi^n + \theta \phi^{n+1} \f$.
   * Generally, only the values 1 and 0.5 are used. The user is not allowed to
   * modify this variable.
   * - 1: first-order
   * - 0.5: second-order\n
   * For the pressure, \ref thetav is always 1. For  the other variables,
   * \ref thetav = 0.5 is used when the  second-order time scheme is activated
   * (\ref optcal::ischtp "ischtp = 2", standard for LES calculations),
   * otherwise \ref thetav = 1.
   *
   * \var blencv
   * Proportion of second-order convective scheme (0 corresponds to an upwind
   * first-order scheme); in case of LES calculation, a second-order scheme is
   * recommended and activated by default (\ref blencv = 1).\n
   * Relevant where \ref iconv = 1.
   *
   * \var blend_st
   * Proportion of second-order convective scheme (0 corresponds to an upwind
   * first-order scheme) after the slope test is activated; in case of LES
   * calculation, a second-order scheme is recommended and activated by
   * default (\ref blend_st = 1).\n
   * Relevant where\ref iconv = 1.
   *
   * \var epsilo
   * Relative precision for the solution of the linear system.  The default is
   * \ref epsilo = \f$ 10^-8 \f$. When there are enough iterations on the
   * reconstruction of the right-hand side of the equation, the value may be
   * increased (by default, in case of second-order in time, with \ref nswrsm =
   * 5 or 10, \ref epsilo is increased to \f$ 10^-5 \f$.
   *
   * \var epsrsm
   * Relative precision on the reconstruction of the right hand-side. The
   * default is \ref epsrsm = \f$ 10^-8 \f$. When there are not enough
   * iterations on the reconstruction of the right-hand side of the equation,
   * the value may be increased (by default, in case of second-order in time,
   * with \ref nswrsm = 5 or 10, \ref epsrsm is increased to \f$ 10^-5 \f$ ).
   *
   * \var epsrgr
   * Relative precision for the iterative gradient reconstruction.
   * (when \ref imrgra = 0).
   *
   * \var climgr
   * For least squares gradients, factor of gradient limitation
   * (high value means little limitation).\n
   * Relevant for all the variables using least-squares gradients for which
   * \ref imligr > CS_GRADIENT_LIMIT_NONE.
   *
   * \var relaxv
   * Relaxation coefficient for the associated variable. This relaxation
   * parameter is only useful for the pressure with the unsteady algorithm (so
   * as to improve the convergence in case of meshes of insufficient quality or
   * of some turbulent models (k-epsilon, v2f, k-omega) and
   * \ref cs_turb_rans_model_t::ikecou "ikecou" = 0; if
   * \ref cs_turb_rans_model_t::ikecou "ikecou" = 1, \ref relaxv is ignored.\n
   * Default values are 0.7 for turbulent variables and 1. for pressure.
   * \ref relaxv also stores the value of the relaxation coefficient when using
   * the steady algorithm, deduced from the value of
   * \ref cs_time_step_options_t::relxst "relxst" (defaulting to
   * \ref relaxv = 1. - relxst).\n
   * Used only for the pressure and for turbulent variables
   * (\f$ k-\epsilon \f$, v2f or \f$ k-\omega \f$ models without coupling) with
   * the unsteady algorithm. Always used with the steady algorithm.
   *
   * \var b_gradient_r
   * Type of boundary gradient reconstruction
   * Same codes as for \ref imrgra
   * (default is 2: least-squares, using extended neighborhood if available)
  */

  int iconv;
  int istat;
  int idircl;
  int ndircl;
  int idiff;
  int idifft;
  int idften;
  int iswdyn;
  int ischcv;
  int ibdtso;
  int isstpc;
  int nswrgr;
  int nswrsm;
  int imvisf;
  int imrgra;
  int imligr;
  int ircflu;
  int iwgrec;
  int icoupl;
  double thetav;   /* TODO: merge with theta */
  double blencv;
  double blend_st;
  double epsilo;
  double epsrsm;
  double epsrgr;
  double climgr;
  double relaxv;

  int  b_gradient_r;

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

  cs_param_bc_type_t          default_bc;
  int                         n_bc_defs;
  cs_xdef_t                 **bc_defs;

  cs_param_bc_enforce_t       default_enforcement;
  cs_real_t                   strong_pena_bc_coeff;
  cs_real_t                   weak_pena_bc_coeff;

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

  int                         n_ic_defs;
  cs_xdef_t                 **ic_defs;

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

  bool                        do_lumping;

  /*!
   * \var time_hodgep
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

  cs_hodge_param_t            time_hodgep;
  cs_property_t              *time_property;
  cs_param_time_scheme_t      time_scheme;
  cs_real_t                   theta;

  /*!
   * @}
   * @name Numerical settings for the diffusion (Laplacian div-grad) term
   * @{
   *
   * \var diffusion_hodgep
   * Set of parameters for the discrete Hodge operator related to the diffusion
   *
   * \var diffusion_property
   * Pointer to the property related to the diffusion term
   */

  cs_hodge_param_t            diffusion_hodgep;
  cs_property_t              *diffusion_property;

  /*!
   * @}
   * @name Numerical settings for the curl-curl term
   * @{
   *
   * \var curlcurl_hodgep
   * Set of parameters for the discrete Hodge operator related to the
   * curl-curl operator
   *
   * \var curlcurl_property
   * Pointer to the property related to the curl-curl term
   */

  cs_hodge_param_t            curlcurl_hodgep;
  cs_property_t              *curlcurl_property;

  /*!
   * @}
   * @name Numerical settings for the grad-div term
   * @{
   *
   * \var graddiv_hodgep
   * Set of parameters for the discrete Hodge operator related to the grad-div
   * operator
   *
   * \var graddiv_property
   * Pointer to the property related to the grad-div term
   */

  cs_hodge_param_t            graddiv_hodgep;
  cs_property_t              *graddiv_property;

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
   * \var adv_strategy
   * Strategy used to handle the advection term (please refer to \ref
   * cs_param_advection_strategy_t)
   *
   * \var adv_extrapol
   * Extrapolation used to estimate the advection field (please refer to \ref
   * cs_param_advection_extrapol_t)
   *
   * \var upwind_portion
   * Value between 0. and 1. (0: centered scheme, 1: pure upwind scheme)
   * Introduce a constant portion of upwinding in a centered scheme
   * Only useful if the advection scheme is set to
   * CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND
   *
   * \var adv_field
   * Pointer to the \ref cs_adv_field_t structure associated to the advection
   * term
   *
   * \var adv_scaling_property
   * May be set to NULL even if the advection term is activated. The value of
   * this property in each cell is a multiplicative coefficient in front of the
   * advection term (boundary terms are also considered) This is useful to treat
   * the thermal module using the variable temperature instead of the enthalpy
   * for instance or in the solidification module.
   */

  cs_param_advection_form_t             adv_formulation;
  cs_param_advection_scheme_t           adv_scheme;
  cs_param_advection_strategy_t         adv_strategy;
  cs_param_advection_extrapol_t         adv_extrapol;
  cs_real_t                             upwind_portion;
  cs_adv_field_t                       *adv_field;
  cs_property_t                        *adv_scaling_property;

  /*!
   * @}
   * @name Numerical settings for the reaction term
   * The contribution of a reaction term to the algebraic system is either at
   * the left-hand and/or right-hand side according to the choice of time
   * scheme
   * @{
   *
   * \var reaction_hodgep
   * Set of parameters for the discrete Hodge operator related to the reaction
   * terms
   *
   * \var n_reaction_terms
   * Number of reaction terms to consider.
   *
   * \var reaction_properties
   * List of properties associated to each reaction term
   */

  cs_hodge_param_t            reaction_hodgep;
  int                         n_reaction_terms;
  cs_property_t             **reaction_properties;

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

  int                         n_source_terms;
  cs_xdef_t                 **source_terms;

  /*!
   * @}
   * @name Definition of the related volume mass injection
   * The contribution of a volume mass injection to the algebraic system
   * is always at right-hand side whatever is the choice of time scheme,
   * and is defined in a manner similar to boundary conditions. For
   * variables whose injection value matches the "ambient" value, no term
   * needs to be added.
   * @{
   *
   * \var n_volume_mass_injections
   * Number of volume injections to consider.
   *
   * \var volume_mass_injections
   * List of definitions of injection values
   */

  int                         n_volume_mass_injections;
  cs_xdef_t                 **volume_mass_injections;

  /*!
   * @}
   * @name Enforcement of values inside the computational domain
   *
   * This is different from the enforcement of boundary conditions but rely on
   * the same algebraic manipulation. Only an algebraic manipulation is used to
   * enforce interior/border DoFs.
   *
   * @{
   *
   * \var n_enforcements
   * Number of enforcements which have been specified
   *
   * \var enforcement_params
   * Array of \ref cs_enforcement_param_t structures storing the settings of
   * each enforcement
   */

  int                         n_enforcements;
  cs_enforcement_param_t    **enforcement_params;

  /*!
   * @}
   * @name Settings related to the resolution of the algebraic system
   * @{
   *
   * \var sles_param
   * Set of parameters to specify how to to solve the algebraic system
   *  - iterative solver
   *  - preconditioner
   *  - tolerance...
   *
   * \var saddle_param
   * Set of parameters to specify how to to solve a saddle-point system.
   * Up to now, only CDO cell-based schemes yield this type of problem. Other
   * saddle-point problems arise when dealing with Stokes or Navier-Stokes
   * equations with a monolithic velocity/pressure coupling. In the latter
   * case, this is directly handled by the structures dedicated to the
   * Navier-Stokes solver.
   *
   * \var incremental_algo_type
   * Type of algorithm used to handle a (non-linear) incremental resolution
   *
   * \var incremental_algo_cvg
   * Set of parameters to drive the convergence of the incremental algorithm
   *
   * \var incremental_relax_factor
   * strictly positive value for the relaxation of the incremental algorithm
   * A recommanded range is between 0.5 and 0.9
   *
   * \var incremental_anderson_param
   * Set of parameters dedicated to the driving of the Anderson acceleration.
   * This is only useful if the incremental_algo_type is set to the Anderson
   * acceleration.
   */

  cs_param_sles_t            *sles_param;
  cs_param_sles_saddle_t     *saddle_param;

  cs_param_nl_algo_t          incremental_algo_type;
  cs_param_sles_cvg_t         incremental_algo_cvg;
  cs_real_t                   incremental_relax_factor;
  cs_iter_algo_param_aac_t    incremental_anderson_param;

  /*! @} */

} cs_equation_param_t;

/*----------------------------------------------------------------------------*/

/*! \enum cs_equation_key_t
 *  \brief List of available keys for setting the parameters of an equation
 *
 * \var CS_EQKEY_ADV_EXTRAPOL
 * Choice in the way to extrapolate the advection field when building the
 * advection operator
 * - "none" (default)
 * - "taylor"
 * - "adams_bashforth"
 * Please refer to \ref cs_param_advection_extrapol_t for more details.
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
 * - "mix_centered_upwind" or "hybrid_centered_upwind"
 *   (cf. \ref CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND)
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
 * \var CS_EQKEY_ADV_STRATEGY
 * Strategy used to handle the advection term
 * - "fully_implicit" or "implicit" (default choice)
 * - "linearized" or "implicit_linear"
 * - "explicit"
 * There is a difference between the two first choices when the advection term
 * induces a non-linearity. In this situation, the first choice implies that a
 * non-linear algorithm has to be used.
 *
 * \var CS_EQKEY_ADV_UPWIND_PORTION
 * Value between 0 and 1 specifying the portion of upwind added to a centered
 * discretization.
 *
 * \var CS_EQKEY_AMG_TYPE
 * Specify which type of algebraic multigrid (AMG) to choose.
 * Available choices are:
 * - "none" --> (default) No predefined AMG solver
 * - "boomer"/"boomer_v"/"bamg" --> Boomer AMG V-cycle multigrid. The HYPRE
                                    library is needed in this case during the
                                    installation process.
 * - "boomer_w"/"bamg_w" --> Boomer AMG W-cycle multigrid. The HYPRE library is
                                    needed in this case during the installation
                                    process.
 * - "gamg"/"gamg_v" --> GAMG V-cycle multigrid. The PETSc library is needed
                         in this case during the installation process.
 * - "gamg_w" --> GAMG W-cycle multigrid from the PETSc library. The PETSc
                  library is needed in this case during the installation
                  process.
 * - "pcmg" --> MG multigrid as preconditioner from the PETSc library
 * - "v_cycle" --> code_saturne's built-in multigrid with a V-cycle strategy

 * - "k_cycle"/ "kamg" --> code_saturne's built-in multigrid with a K-cycle
 *                         strategy (cf. Notay's articles on that subject).
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
 * - "de_rham" or "derham" (default): Evaluation at vertices for potentials,
 *   integral along a line for circulations, integral across the normal
 *   component of a face for fluxes and integral inside a cell for
 *   densities. In case of "potential", no quadrature rule is applied.
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
 *                 consistent for all meshes. Require an "orthogonal" (or
 *                 admissible) mesh; Warning: do not use this algorithm on
 *                 general meshes.
 * - "bubble" --> Reconstruction op. relying on a stabilization which is built
 *                with a bubble function associated to edges

 * - "cost" or "ocs" --> Orthogonal decomposition of the Consistent and
 *                       Stabilized part. (cost = COnsitent + STabilized) The
 *                       stabilization part of the reconstruction op. is a
 *                       piecewise constant polynom on a partition of each cell
 * - "wbs" --> is robust and accurate but is limited to the reconstruction of
 *             potential-like degrees of freedom and needs a correct
 *             computation of the cell barycenter. Make sense only in the case
 *             of a CDOVB or CDOVCB scheme. In the case of tetrahedral meshes,
 *             the stiffness matrix is the same as the one encountered in FE P1
 *             schemes
 *
 * \var CS_EQKEY_HODGE_DIFF_COEF
 * This key is only useful if CS_EQKEY_HODGE_{TIME, DIFF, REAC}_ALGO is set to
 * "cost".
 * keyval is either a name or a value:
 * - "dga" corresponds to the value \f$ 1./3. \f$
 * - "sushi" corresponds to the value \f$1./\sqrt(3.)\f$
 * - "gcr"  corresponds to the value \f$1\f$.
 * - "frac23" or "2/3" corresponds to the value \f$2./3.\f$.
 * - use an expression like "1.5", "9" for a user-defined value
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
 * equation. Avalaible choices are:\n
 * - "amg"           --> Algebraic MultiGrid iterative solver.
 *                       Good choice for a scalable solver related to symmetric
 *                       positive definite system.
 * - "jacobi","diag"
 *    or "diagonal"  --> Jacobi algorithm: simpliest iterative solver
 * - "gauss_seidel"
 *    or "gs"        --> Gauss-Seidel algorithm
 * - "cg"            --> conjuguate gradient algorithm
 * - "cr3"           --> a 3-layer conjugate residual solver (when "cs" is
 *                       chosen as the solver family)
 * - "fcg"           --> flexible version of the conjuguate gradient
 *                       algorithm used when the preconditioner can change
 *                       iteration by iteration
 * - "bicg"          --> Bi-CG algorithm (for non-symmetric linear systems)
 * - "bicgstab2"     --> BiCG-Stab2 algorithm (for non-symmetric linear
 *                       systems)
 * - "gcr"           --> robust and flexible iterative solver. Not the best
 *                       choice if the system is easy to solve
 * - "gmres"         --> robust iterative solver. Not the best choice if the
 *                       system is easy to solve
 * - "fgmres"        --> Flexible gmres (only with PETSc installation up to
 *                       now). An evolutive preconditioner can be used with
 *                       this solver. This is a very robust iterative solver.
 *                       Not the best choice if the system is easy to solve
 * - "minres"        --> Solver of choice for symmetric indefinite systems
 * - "mumps"         --> Direct solver (MUMPS or PETSc)
 *                       By default, it performs a LU factorization on general
 *                       matrices. More options are available using the function
 *                       \ref cs_param_sles_mumps and the function
 *                       \ref cs_param_sles_mumps_advanced
 * - "sym_gauss_seidel" or "sgs" --> Symmetric Gauss-Seidel algorithm
 * - "user"          --> User-defined iterative solver (rely on the function
 *                       cs_user_sles_it_solver())
 * - "none"          --> No solver.
 *
 * \var CS_EQKEY_ITSOL_ATOL

 * Absolute tolerance factor for stopping the iterative process during the
 * iterative resolution of a linear system related to an equation. Most
 * iterative solver are not using this tolerance (the relative tolerance is
 * always used). PETSc solvers use this information for instance.Please refer
 * to \ref CS_EQKEY_ITSOL_RTOL\n
 * - Example: "1e-14"
 *
 * \var CS_EQKEY_ITSOL_DTOL
 * Divergence tolerance factor for stopping the iterative process during the
 * iterative resolution of a linear system related to an equation. Most
 * iterative solver are not using this tolerance (the relative tolerance is
 * always used). PETSc solvers use this information for instance.Please refer
 * to \ref CS_EQKEY_ITSOL_RTOL\n
 * - Example: "1e3"
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
 * "filtered_rhs" or "filtered_rhs"
 *
 * \var CS_EQKEY_ITSOL_RESTART
 * Maximum number of iterations before restarting a Krylov solver
 * Only useful with GMRES, flexible GMRES or GCR solvers.
 * This value has an impact on the memory footprint since one has to store a
 * buffer of double with a size equal to restart*sizeof(solution array)
 * - Example: "20"
 *
 * \var CS_EQKEY_ITSOL_RTOL
 * \var CS_EQKEY_ITSOL_EPS (deprecated)
 * Relative tolerance factor for stopping the iterative process during the
 * iterative resolution of a linear system related to an equation.\n
 * - Example: "1e-10"
 *
 * \var CS_EQKEY_PRECOND
 * Specify the preconditioner associated to an iterative solver. Be careful
 * some options are only available with a given solver class. Be sure that your
 * installation has been installed with the appropriate library.\n
 * Available choices are:
 * - "none": no preconditioner is used
 * - "jacobi" or "diag": diagonal preconditoner
 * - "block_jacobi"/"bjacobi": Block Jacobi with ILU(0) in each block. A block
 *                             is the matrix associated to a rank. (Only with
 *                             PETSc)
 * - "bjacobi_sgs"/"bjacobi_ssor": Block Jacobi with Symm. Gauss-Seidel in each
 *                             block. A block is the matrix associated to a
 *                             rank. (Only with PETSc)
 * - "poly1": Neumann polynomial of order 1 (only with code_saturne)
 * - "poly2": Neumann polynomial of order 2 (only with code_saturne)
 * - "ssor": symmetric successive over-relaxation (only with PETSC)
 * - "ilu0": incomplete LU factorization (only with PETSc)
 * - "icc0": incomplete Cholesky factorization (for symmetric matrices and
 *           only with PETSc)
 * - "lu": LU factorization (only with PETSc). It may use MUMPS if PETSc is
 *         built with MUMPS
 * - "amg": algebraic multigrid technique (see \ref CS_EQKEY_AMG_TYPE for
            additional settings)
 * - "amg_block"/"block_amg: algebraic multigrid by block (useful for
 *                vector-valued equations). By default, a diagonal block
 *                preconditioning is used if nothing else is set.

 * - "mumps": Direct solver (only MUMPS)
 *            By default, it performs a LU factorization on general matrices.
 *            More options are available using the functions
 *            \ref cs_param_sles_mumps and \ref cs_param_sles_mumps_advanced
 *
 * \var CS_EQKEY_PRECOND_BLOCK_TYPE
 * Specify the type of block preconditioner associated to a preconditioner.
 * When "full" is specified. That means that the smallest possible blocks are
 * considered. Be careful: Most of these options are available only with
 * PETSc/HYPRE.\n
 * Available choices are:
 * - "none": no block preconditioner (default choice)
 * - "diag": diagonal (or additive) block preconditoner
 * - "lower": lower triangular (multiplicative) block preconditioner
 * - "upper": upper triangular (multiplicative) block preconditioner
 * - "symm": symmetric Gauss-Seidel block preconditioner
 * - "full_diag": full diagonal (or additive) block preconditoner.
 * - "full_lower": full lower triangular (multiplicative) block preconditioner
 * - "full_upper": upper triangular (multiplicative) block preconditioner
 * - "full_symm": symmetric Gauss-Seidel block preconditioner
 *
 * \var CS_EQKEY_SADDLE_PRECOND
 * Preconditioner used to solve a saddle-point system.\n
 * Please refer to \ref cs_param_saddle_precond_t for more details.\n
 * Available choices are:\n
 * - "none" No preconditioner is used (default)
 * - "diag"
 * - "lower"
 * - "upper"
 *
 * \var CS_EQKEY_SADDLE_RTOL
 * Relative tolerance under which the iterative process stops
 * - Example: "1e-7"
 *
 * \var CS_EQKEY_SADDLE_SOLVER
 * Strategy/solver used to solve a saddle-point system.\n
 * Available choices are:
 * - "none" No solver (default --> no saddle-point system)
 * - "gcr"  GCR iterative solver
 * - "mumps" MUMPS sparse direct solver (external library)
 * - "minres" MINRES iterative solver
 *
 * \var CS_EQKEY_SADDLE_VERBOSITY
 * Level of details displayed for the resolution of a saddle-point system
 * - Examples: "0", "1", "2" or higher
 *
 * \var CS_EQKEY_SCHUR_APPROX
 * Choice of apprixmation of the Schur complement. Only useful when a
 * saddle-point problem has to be solved. Some solvers or choices of
 * preconditioner/solver may not need an approximation of the schur complement.
 * Please refer to \ref cs_param_schur_approx_t for more details.\n
 * Available choices are:
 * "none" no Schur approximation --> the identity is used
 * "diag_inv" use the inverse of the (1,1)-block diagonal
 * "lumped_inv" use the inverse of the (1,1)-block after lumping
 * "mass" use a scaled mass matrix for the space related to the (2,2)-block
 *
 * \var CS_EQKEY_SLES_VERBOSITY
 * Level of details written by the code for the resolution of the linear system
 * - Examples: "0", "1", "2" or higher
 *
 * \var CS_EQKEY_SPACE_SCHEME
 * Set the space discretization scheme. Available choices are:
 * - "cdo_vb"  or "cdovb" for CDO vertex-based scheme
 * - "cdo_vcb" or "cdovcb" for CDO vertex+cell-based scheme
 * - "cdo_fb"  or "cdofb" for CDO face-based scheme
 * - "cdo_eb"  or "cdoeb" for CDO edge-based scheme
 * - "cdo_cb"  or "cdocb" for CDO cell-based scheme
 * - "hho_p1" for HHO schemes with \f$\mathbb{P}_1\f$ polynomial approximation
 * - "hho_p2" for HHO schemes with \f$\mathbb{P}_2\f$ polynomial approximation
 *
 * \var CS_EQKEY_TIME_SCHEME
 * Set the scheme for the temporal discretization. Available choices are:
 * - "euler_implicit" or "forward_euler" (1st order)
 * - "euler_explicit" or "backward_euler" (1st order, conditional stability)
 * - "theta_scheme": Time scheme encompassing several other schemes according
 *   to the value of theta. One recovers for instance "euler_implicit" with
 *   theta equal to "1", "euler_explicit" with "0" or "crank_nicolson" with
 *   theta equal to "0.5"
 * - "crank_nicolson": Shortcut to set a theta scheme with theta equal to
 *    0.5. This time discretization is second_order in time accurate.
 * - "bdf2": 2nd order, implicit time scheme
 *
 * \var CS_EQKEY_TIME_THETA
 * Set a value betwwen 0 and 1 for the theta parameter when a "theta_scheme" is
 * set for the time discretization. Only useful if CS_EQKEY_TIME_SCHEME is set
 * to "theta_scheme".
 * - Example: "0.75"
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

  CS_EQKEY_ADV_EXTRAPOL,
  CS_EQKEY_ADV_FORMULATION,
  CS_EQKEY_ADV_SCHEME,
  CS_EQKEY_ADV_STRATEGY,
  CS_EQKEY_ADV_UPWIND_PORTION,
  CS_EQKEY_AMG_TYPE,
  CS_EQKEY_BC_ENFORCEMENT,
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
  CS_EQKEY_ITSOL_ATOL,
  CS_EQKEY_ITSOL_DTOL,
  CS_EQKEY_ITSOL_EPS,
  CS_EQKEY_ITSOL_MAX_ITER,
  CS_EQKEY_ITSOL_RESNORM_TYPE,
  CS_EQKEY_ITSOL_RESTART,
  CS_EQKEY_ITSOL_RTOL,
  CS_EQKEY_OMP_ASSEMBLY_STRATEGY,
  CS_EQKEY_PRECOND,
  CS_EQKEY_PRECOND_BLOCK_TYPE,
  CS_EQKEY_SADDLE_PRECOND,
  CS_EQKEY_SADDLE_RTOL,
  CS_EQKEY_SADDLE_SOLVER,
  CS_EQKEY_SADDLE_VERBOSITY,
  CS_EQKEY_SCHUR_APPROX,
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
 * \brief  Ask if the parameters of the equation induces an implicit treatment
 *         of the advection term
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_implicit_advection(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_CONVECTION) {
    if (eqp->adv_strategy == CS_PARAM_ADVECTION_IMPLICIT_FULL ||
        eqp->adv_strategy == CS_PARAM_ADVECTION_IMPLICIT_LINEARIZED)
      return true;
    else
      return false;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameters of the equation has activated a user hook
 *         to get a fine tuning of the cellwise system building
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_equation_param_has_user_hook(const cs_equation_param_t     *eqp)
{
  assert(eqp != NULL);
  if (eqp->flag & CS_EQUATION_BUILD_HOOK)
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
cs_equation_param_create(const char            *name,
                         cs_equation_type_t     type,
                         int                    dim,
                         cs_param_bc_type_t     default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a \ref cs_equation_param_t structure
 *
 * \deprecated Renamed to\ref cs_equation_param_create.
 *
 * \param[in] name          name of the equation
 * \param[in] type          type of equation
 * \param[in] dim           dim of the variable associated to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return a pointer to a new allocated \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

inline static cs_equation_param_t *
cs_equation_create_param(const char            *name,
                         cs_equation_type_t     type,
                         int                    dim,
                         cs_param_bc_type_t     default_bc)
{
  return cs_equation_param_create(name, type, dim, default_bc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the settings from one \ref cs_equation_param_t structure to
 *         another one. The name is not copied.
 *
 * \param[in]      ref          pointer to the reference cs_equation_param_t
 * \param[in, out] dst          pointer to the cs_equation_param_t to update
 * \param[in]      copy_fld_id  copy also the field id or not
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_copy_from(const cs_equation_param_t   *ref,
                            cs_equation_param_t         *dst,
                            bool                         copy_fld_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy only the part dedicated to the boundary conditions and the DoF
 *        (degrees of freedom) enforcement from one \ref cs_equation_param_t
 *        structure to another one.
 *
 * \param[in]      ref       pointer to the reference \ref cs_equation_param_t
 * \param[in, out] dst       pointer to the \ref cs_equation_param_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_copy_bc(const cs_equation_param_t   *ref,
                          cs_equation_param_t         *dst);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the settings from one \ref cs_equation_param_t structure to
 *         another one. The name is not copied.
 *
 * \deprecated Renamed to\ref cs_equation_param_copy_from.
 *
 * \param[in]      ref       pointer to the reference \ref cs_equation_param_t
 * \param[in, out] dst       pointer to the \ref cs_equation_param_t to update
 * \param[in]      copy_fid  copy also the field id or not
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_equation_copy_param_from(const cs_equation_param_t   *ref,
                            cs_equation_param_t         *dst,
                            bool                         copy_fid)
{
  cs_equation_param_copy_from(ref, dst, copy_fid);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the contents of a \ref cs_equation_param_t
 *
 * The cs_equation_param_t structure itself is not freed, but all its
 * sub-structures are freed.
 *
 * This is useful for equation parameters which are accessed through
 * field keywords.
 *
 * \param[in, out]  eqp  pointer to a \ref cs_equation_param_t
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_clear(cs_equation_param_t   *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the contents of a \ref cs_equation_param_t
 *
 * \deprecated Renamed to\ref cs_equation_param_clear.
 *
 * The cs_equation_param_t structure itself is not freed, but all its
 * sub-structures are freed.
 *
 * This is useful for equation parameters which are accessed through
 * field keywords.
 *
 * \param[in, out]  eqp  pointer to a \ref cs_equation_param_t
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_equation_clear_param(cs_equation_param_t   *eqp)
{
  cs_equation_param_clear(eqp);
}

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
cs_equation_param_free(cs_equation_param_t     *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_equation_param_t
 *
 * \deprecated Renamed to\ref cs_equation_param_free.
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

inline static cs_equation_param_t *
cs_equation_free_param(cs_equation_param_t     *eqp)
{
  return cs_equation_param_free(eqp);
}

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
cs_equation_param_set(cs_equation_param_t   *eqp,
                      cs_equation_key_t      key,
                      const char            *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a \ref cs_equation_param_t
 *         structure
 *
 * \deprecated Renamed to\ref cs_equation_param_set.
 *
 * \param[in, out]  eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_equation_set_param(cs_equation_param_t   *eqp,
                      cs_equation_key_t      key,
                      const char            *keyval)
{
  cs_equation_param_set(eqp, key, keyval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle the SLES solver
 *        associated to this set of equation parameters
 *
 * \param[in] eqp      pointer to a \ref cs_equation_param_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_equation_param_get_sles_param(cs_equation_param_t   *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Settings are related to this equation.
 *
 * \param[in, out] eqp         pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_sles(cs_equation_param_t      *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the given quadrature rule to all existing definitions under the
 *        cs_equation_param_t structure. To get a more detailed control of the
 *        quadrature rule, please consider the function
 *        \ref cs_xdef_set_quadrature
 *
 * \param[in, out] eqp     pointer to a \ref cs_equation_param_t structure
 * \param[in]      qtype   type of quadrature to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_quadrature_to_all(cs_equation_param_t   *eqp,
                                        cs_quadrature_type_t   qtype);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lock settings to prevent from unwanted modifications.
 *
 * \param[in, out] eqp   pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_lock_settings(cs_equation_param_t   *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Unlock settings. Be sure that is really wanted (inconsistency between
 *        the setup logging and what is used may appear)
 *
 * \param[in, out] eqp   pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_unlock_settings(cs_equation_param_t   *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief At this stage, the numerical settings should not be modified anymore
 *        by the user. One makes a last set of modifications if needed to
 *        ensure a consistent numerical settings.
 *
 * \param[in, out] eqp      pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_ensure_consistent_settings(cs_equation_param_t   *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the detail of a \ref cs_equation_param_t structure
 *
 * \param[in]  eqp      pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_log(const cs_equation_param_t  *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameter settings of the equation has requested the
 *         treatment of Robin boundary conditions
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_param_has_robin_bc(const cs_equation_param_t     *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this
 *         equation. This definition applies to a volume zone.
 *         By default, the unknown is set to zero everywhere.
 *         Here a constant value is set to all the unknows belonging to the
 *         given zone with name z_name
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
 * \brief  Define the initial condition for the unknown related to this
 *         equation. This definition applies to a volume zone.
 *         By default, the unknown is set to zero everywhere.
 *         Here the value set to each unknown belonging to the given zone with
 *         name z_name is such that the integral over the cells of the zone
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
 *         equation. This definition applies to a volume zone.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value for each unknown associated to the zone with
 *         name z_name is set according to an analytical function
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
 * \brief  Define the initial condition for the unknown related to this
 *         equation. This definition applies to a volume zone.
 *         By default, the unknown is set to zero everywhere.
 *         Case of a definition by a DoF function.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  where information is computed
 * \param[in]      func      pointer to a DoF function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_dof_func(cs_equation_param_t    *eqp,
                               const char             *z_name,
                               cs_flag_t               loc_flag,
                               cs_dof_func_t          *func,
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
 * \param[in, out] eqp           pointer to a cs_equation_param_t structure
 * \param[in]      bc_type       type of boundary condition to add
 * \param[in]      z_name        name of the related boundary zone
 * \param[in]      loc           information to know where are located values
 * \param[in]      array         pointer to an array
 * \param[in]      is_owner      transfer the lifecycle to the cs_xdef_t struct.
 *                               (true or false)
 * \param[in]      full_length   if true, size of "array" should be allocated
 *                               to the total numbers of entities related to the
 *                               given location. If false, a new list is
 *                               allocated and filled with the related subset
 *                               indirection.
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
                            bool                        full_length);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       field     pointer to a cs_field_t structure
 *
 * \return a pointer to the new allocated \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_field(cs_equation_param_t        *eqp,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_field_t                 *field);

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
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given equation param structure
 *        ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *        Definition relying on a \ref cs_time_func_t function pointer
 *
 * \param[in, out] eqp      pointer to a cs_equation_param_t structure
 * \param[in]      bc_type  type of boundary condition to add
 * \param[in]      z_name   name of the associated zone (if NULL or "" if
 *                          all cells are considered)
 * \param[in]      t_func   pointer to an analytic function defining the value
 * \param[in]      input    NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_time_func(cs_equation_param_t        *eqp,
                                const cs_param_bc_type_t    bc_type,
                                const char                 *z_name,
                                cs_time_func_t             *t_func,
                                void                       *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given cs_equation_param_t structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  location where values are computed
 * \param[in]      func      pointer to cs_dof_func_t function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_dof_func(cs_equation_param_t        *eqp,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *z_name,
                               cs_flag_t                   loc_flag,
                               cs_dof_func_t              *func,
                               void                       *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return pointer to existing boundary condition definition structure
 *         for the given equation param structure and zone.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 *
 * \return a pointer to the \ref cs_xdef_t structure if present, or NULL
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_find_bc(cs_equation_param_t   *eqp,
                    const char            *z_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Remove boundary condition from the given equation param structure
 *         for a given zone.
 *
 * If no matching boundary condition is found, the function returns
 * silently.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_remove_bc(cs_equation_param_t   *eqp,
                      const char            *z_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Remove initial condition from the given equation param structure
 *         for a given zone.
 *
 * If no matching boundary condition is found, the function returns
 * silently.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_remove_ic(cs_equation_param_t   *eqp,
                      const char            *z_name);

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
 * \param[in]      inversion  1: true, false otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_curlcurl(cs_equation_param_t   *eqp,
                         cs_property_t         *property,
                         int                    inversion);

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
 * \brief Associate a new term related to the time derivative operator for
 *        the equation associated to the given \ref cs_equation_param_t
 *        structure
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
 * \brief Associate a new term related to the advection operator for the
 *        equation associated to the given \ref cs_equation_param_t structure
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
 * \brief Associate a scaling property to the advection
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_advection_scaling_property(cs_equation_param_t   *eqp,
                                           cs_property_t         *property);

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
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      z_name       name of the associated zone (if NULL or "" if
 *                              all cells are considered)
 * \param[in]      loc          information to know where are located values
 * \param[in]      array        pointer to an array
 * \param[in]      is_owner     transfer the lifecycle to the cs_xdef_t struct.
 *                              (true or false)
 * \param[in]      full_length  if true, the size of "array" should be allocated
 *                              to the total numbers of entities related to the
 *                              given location. If false, a new list is
 *                              allocated and filled with the related subset
 *                              indirection.
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
                                     bool                    full_length);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using a constant value.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_value(cs_equation_param_t  *eqp,
                                               const char           *z_name,
                                               double               *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using a constant quantity
 *         distributed over the associated zone's volume.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      quantity  pointer to quantity to distribute over the zone
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_qov(cs_equation_param_t  *eqp,
                                             const char           *z_name,
                                             double               *quantity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using an analytical function.
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
cs_equation_add_volume_mass_injection_by_analytic(cs_equation_param_t   *eqp,
                                                  const char            *z_name,
                                                  cs_analytic_func_t    *func,
                                                  void                  *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using a DoF function.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  where information is computed
 * \param[in]      func      pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_dof_func(cs_equation_param_t  *eqp,
                                                  const char           *z_name,
                                                  cs_flag_t             loc_flag,
                                                  cs_dof_func_t        *func,
                                                  void                 *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value of degrees of freedom located at
 *         the mesh vertices.
 *         The spatial discretization scheme for the given equation has to be
 *         CDO vertex-based or CDO vertex+cell-based schemes.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      n_vertices   number of vertices to enforce
 * \param[in]      vertex_ids   list of vertices
 * \param[in]      ref_value    default values or ignored (may be NULL)
 * \param[in]      vtx_values   list of associated values, ignored if NULL
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_vertex_dof_enforcement(cs_equation_param_t    *eqp,
                                       cs_lnum_t               n_vertices,
                                       const cs_lnum_t         vertex_ids[],
                                       const cs_real_t         ref_value[],
                                       const cs_real_t         vtx_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value of degrees of freedom located at
 *         the mesh edges.
 *         The spatial discretization scheme for the given equation has to be
 *         CDO edge-based schemes.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp           pointer to a cs_equation_param_t structure
 * \param[in]      n_edges       number of edges to enforce
 * \param[in]      edge_ids      list of edges
 * \param[in]      ref_value     default values or ignored (may be NULL)
 * \param[in]      edge_values   list of associated values, ignored if NULL
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_edge_dof_enforcement(cs_equation_param_t    *eqp,
                                     cs_lnum_t               n_edges,
                                     const cs_lnum_t         edge_ids[],
                                     const cs_real_t         ref_value[],
                                     const cs_real_t         edge_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value of degrees of freedom located at
 *         the mesh faces.
 *         The spatial discretization scheme for the given equation has to be
 *         CDO face-based schemes.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp           pointer to a cs_equation_param_t structure
 * \param[in]      n_faces       number of faces to enforce
 * \param[in]      face_ids      list of faces
 * \param[in]      ref_value     default values or ignored (may be NULL)
 * \param[in]      face_values   list of associated values, ignored if NULL
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_face_dof_enforcement(cs_equation_param_t    *eqp,
                                     cs_lnum_t               n_faces,
                                     const cs_lnum_t         face_ids[],
                                     const cs_real_t         ref_value[],
                                     const cs_real_t         face_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value related to the degrees of freedom
 *         associated to the list of selected cells.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      n_cells      number of selected cells
 * \param[in]      cell_ids     list of cell ids
 * \param[in]      ref_value    ignored if NULL
 * \param[in]      cell_values  list of associated values, ignored if NULL
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_cell_enforcement(cs_equation_param_t   *eqp,
                                 cs_lnum_t              n_cells,
                                 const cs_lnum_t        elt_ids[],
                                 const cs_real_t        ref_value[],
                                 const cs_real_t        cell_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new enforcement if enforcement_id does not exist or replace it
 *         otherwise. Enforcement of the value related to the degrees of freedom
 *         associated to the list of selected cells.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp             pointer to a cs_equation_param_t structure
 * \param[in]      enforcement_id  id of the enforcement to handle
 * \param[in]      n_cells         number of selected cells
 * \param[in]      cell_ids        list of cell ids
 * \param[in]      ref_value       ignored if NULL
 * \param[in]      cell_values     list of associated values, ignored if NULL
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_or_replace_cell_enforcement(cs_equation_param_t *eqp,
                                            int                  enforcement_id,
                                            cs_lnum_t            n_cells,
                                            const cs_lnum_t      cell_ids[],
                                            const cs_real_t      ref_value[],
                                            const cs_real_t      cell_values[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_PARAM_H__ */
