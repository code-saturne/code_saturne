#ifndef __CS_PARAM_H__
#define __CS_PARAM_H__

/*============================================================================
 * Manage the definition/setting of a computation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/* Tag to build parameter flag */
#define CS_PARAM_FLAG_UNIFORM  (1 <<  0) //    1: uniform (in space)
#define CS_PARAM_FLAG_CELLWISE (1 <<  1) //    2: cellwise uniform
#define CS_PARAM_FLAG_UNSTEADY (1 <<  2) //    4: unsteady
#define CS_PARAM_FLAG_VERTEX   (1 <<  3) //    8: on vertices
#define CS_PARAM_FLAG_EDGE     (1 <<  4) //   16: on edges
#define CS_PARAM_FLAG_FACE     (1 <<  5) //   32: on faces
#define CS_PARAM_FLAG_CELL     (1 <<  6) //   64: on cells
#define CS_PARAM_FLAG_PRIMAL   (1 <<  7) //  128: on primal mesh
#define CS_PARAM_FLAG_DUAL     (1 <<  8) //  256: on dual mesh
#define CS_PARAM_FLAG_BORDER   (1 <<  9) //  512: scalar-valued
#define CS_PARAM_FLAG_SCAL     (1 << 10) // 1024: scalar-valued
#define CS_PARAM_FLAG_VECT     (1 << 11) // 2048: vector-valued
#define CS_PARAM_FLAG_TENS     (1 << 12) // 4096: tensor-valued
#define CS_PARAM_FLAG_BY_CELL  (1 << 13) // 8192: by cell (c2e, c2f, c2v)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* User-defined function */
typedef void (cs_user_func_t) (const void         *input1,
                               const void         *input2,
                               cs_real_t           tcur,
                               const cs_real_3_t   xyz,
                               cs_get_t           *output);

typedef union {

  cs_get_t               get;       // definition by value
  cs_analytic_func_t    *analytic;  // definition by an analytic function
  cs_timestep_func_t    *time_func; // definition of the time step by a function
  cs_user_func_t        *user_func; // definition by an user-defined function
  cs_onevar_law_func_t  *law1_func; // definition by law depending one variable

} cs_def_t;

typedef enum {

  CS_PARAM_DEF_BY_ANALYTIC_FUNCTION,
  CS_PARAM_DEF_BY_ARRAY,
  CS_PARAM_DEF_BY_FIELD,
  CS_PARAM_DEF_BY_ONEVAR_LAW,
  CS_PARAM_DEF_BY_TIME_FUNCTION,
  CS_PARAM_DEF_BY_USER_FUNCTION,
  CS_PARAM_DEF_BY_VALUE,
  CS_PARAM_N_DEF_TYPES

} cs_param_def_type_t;

/* Dimension of the variable to deal with */
typedef enum {

  CS_PARAM_VAR_SCAL,    // scalar variable (dim = 1)
  CS_PARAM_VAR_VECT,    // vector variable (dim = 3)
  CS_PARAM_VAR_SYMTENS, // symmetric tensor variable (dim = 6)
  CS_PARAM_VAR_TENS,    // tensor variable (dim = 9)
  CS_PARAM_N_VAR_TYPES

} cs_param_var_type_t;

/* DISCRETE HODGE OPERATORS */
/* ======================== */

typedef enum {

  CS_PARAM_HODGE_TYPE_VPCD, // from primal vertices to dual cells
  CS_PARAM_HODGE_TYPE_EPFD, // from primal edges to dual faces
  CS_PARAM_HODGE_TYPE_FPED, // from primal faces to dual edges
  CS_PARAM_HODGE_TYPE_EDFP, // from dual edges to primal faces
  CS_PARAM_HODGE_TYPE_CPVD, // from primal cells to dual vertices
  CS_PARAM_N_HODGE_TYPES

} cs_param_hodge_type_t;

typedef enum {

  CS_PARAM_HODGE_ALGO_VORONOI, // Under othogonality condition gives a diag. op.
  CS_PARAM_HODGE_ALGO_WBS,     // WBS: Whitney Barycentric Subdivision
  CS_PARAM_HODGE_ALGO_COST,    // COST: COnsistency & STabilization splitting
  CS_PARAM_N_HODGE_ALGOS

} cs_param_hodge_algo_t;

typedef struct {

  bool                    inv_pty; /* Definition based on the property
                                      or its inverse */

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

/* Parameters related to time discretization */
typedef struct {

  cs_time_scheme_t   scheme;     // numerical scheme
  cs_real_t          theta;      // used in theta-scheme
  bool               do_lumping; // perform mass lumping ?

  /* Initial conditions */
  cs_param_def_type_t   ic_def_type;  // type of definition
  cs_def_t              ic_def;       // definition

} cs_param_time_t;

/* ADVECTION OPERATOR PARAMETRIZATION */
/* ================================== */

typedef enum {

  CS_PARAM_ADVECTION_FORM_CONSERV,
  CS_PARAM_ADVECTION_FORM_NONCONS,
  CS_PARAM_N_ADVECTION_FORMULATIONS

} cs_param_advection_form_t;

/* Type of weight functions considered */
typedef enum {

  CS_PARAM_ADVECTION_WEIGHT_ALGO_CENTERED,
  CS_PARAM_ADVECTION_WEIGHT_ALGO_UPWIND,
  CS_PARAM_ADVECTION_WEIGHT_ALGO_SAMARSKII,
  CS_PARAM_ADVECTION_WEIGHT_ALGO_SG,
  CS_PARAM_ADVECTION_WEIGHT_ALGO_D10G5,
  CS_PARAM_N_ADVECTION_WEIGHT_ALGOS

} cs_param_advection_weight_algo_t;

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

  cs_param_advection_form_t         formulation; // conservative or not
  cs_param_advection_weight_algo_t  weight_algo;
  cs_param_advection_weight_t       weight_criterion;
  cs_quadra_type_t                  quad_type; // barycentric, higher, highest

} cs_param_advection_t;

/* REACTION TERM PARAMETRIZATION */
/* ============================= */

typedef enum {

  CS_PARAM_REACTION_TYPE_LINEAR,
  CS_PARAM_N_REACTION_TYPES

} cs_param_reaction_type_t;

typedef struct {

  char                      *name;
  cs_param_reaction_type_t   type;

  cs_param_hodge_t           hodge;
  bool                       do_lumping;

} cs_param_reaction_t;

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

/* Definition of the value of a boundary condition (BC)
   One or two values may be needed according to the type of BC.
   That's why two coefficients are used in the definition.
   For Dirichlet or Neumann BC for instance, only one coefficient is used.
*/

typedef struct {

  int                    loc_id;   // Id related to the list of border faces

  cs_param_bc_type_t     bc_type;  // type of mathematical BC
  cs_param_var_type_t    var_type; // type of variable
  cs_param_def_type_t    def_type; // Type of definition for a and b


  /* Access to the value related to the first coefficient and possibly
     the second one */
  cs_def_t               def_coef1;
  cs_def_t               def_coef2;

} cs_param_bc_def_t;

typedef struct {

  cs_param_bc_type_t     default_bc;   // BC used by default
  cs_param_bc_enforce_t  enforcement;  // type of boundary enforcement
  cs_quadra_type_t       quad_type;    // barycentric, higher, highest...
  bool                   use_subdiv;   /* subdivide or not into tetrahedra for
                                          evaluating BC values */

  int                    n_defs;
  cs_param_bc_def_t     *defs;

} cs_param_bc_t;

/* SOURCE TERMS */
/* ============ */

/* Types of source terms */
typedef enum {

  CS_PARAM_SOURCE_TERM_USER,     // user-defined
  CS_PARAM_SOURCE_TERM_MASS,     // specific treatment
  CS_PARAM_SOURCE_TERM_HEADLOSS, // specific treatment
  CS_PARAM_N_SOURCE_TERM_TYPES

} cs_param_source_term_type_t;

typedef struct {

  char  *restrict name;   /* short description of the source term */

  int             ml_id;  /* id of the related mesh location structure */
  int             post;   /* -1: no post, 0: at the beginning,
                              n: at each 'n' iterations */

  /* Specification related to the way of computing the source term */
  cs_param_source_term_type_t  type;       // mass, head loss...
  cs_param_var_type_t          var_type;   // scalar, vector...
  cs_param_def_type_t          def_type;   // by value, by function...
  cs_quadra_type_t             quad_type;  // barycentric, higher, highest
  bool                         use_subdiv; // use a subdivision into tetrahedra

  /* Values if one needs an implicit and explicit part
     implicit part comes first */
  cs_def_t                     def;

} cs_param_source_term_t;

/* ITERATIVE SOLVERS */
/* ================= */

typedef enum {

  CS_PARAM_PRECOND_DIAG,    // Diagonal preconditioning (also called Jacobi)
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

  CS_PARAM_ITSOL_CG,        // Conjuguate Gradient
  CS_PARAM_ITSOL_BICG,      // Bi-Conjuguate gradient
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
 * \brief  Set a cs_def_t structure
 *
 * \param[in]      def_type   type of definition (by value, function...)
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      val        value to set
 * \param[in, out] def        pointer to cs_def_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_set_def(cs_param_def_type_t      def_type,
                 cs_param_var_type_t      var_type,
                 const void              *val,
                 cs_def_t                *def);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_param_bc_t structure
 *
 * \param[in]  default_bc     default boundary condition
 *
 * \return a pointer to the new structure (free with cs_param_eq_t)
 */
/*----------------------------------------------------------------------------*/

cs_param_bc_t *
cs_param_bc_create(cs_param_bc_type_t  default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_param_bc_def_t structure
 *
 * \param[in, out] bcpd       pointer to cs_param_bc_def_t struct. to set
 * \param[in]      loc_id     id related to a cs_mesh_location_t
 * \param[in]      bc_type    generic type of admissible boundary conditions
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      def_type   by value, function...
 * \param[in]      coef1      access to the value of the first coef
 * \param[in]      coef2      access to the value of the second coef (optional)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_bc_def_set(cs_param_bc_def_t      *bcpd,
                    int                     loc_id,
                    cs_param_bc_type_t      bc_type,
                    cs_param_var_type_t     var_type,
                    cs_param_def_type_t     def_type,
                    const void             *coef1,
                    const void             *coef2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new reaction term. The structure related to this reaction
 *         term has already been allocated among the list of reaction terms
 *         associated to an equation
 *
 * \param[in, out] rp         pointer to cs_param_reaction_t structure
 * \param[in]      r_name     name of the reaction term
 * \param[in]      h_type     type of discrete Hodge op. associated to this term
 * \param[in]      h_algo     algorithm used to build the discrete Hodge op.
 * \param[in]      r_type     type of reaction term
 */
/*----------------------------------------------------------------------------*/

void
cs_param_reaction_add(cs_param_reaction_t          *rp,
                      const char                   *r_name,
                      cs_param_hodge_type_t         h_type,
                      cs_param_hodge_algo_t         h_algo,
                      cs_param_source_term_type_t   r_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a source term. This source term is added to the list of
 *         source terms associated to an equation
 *
 * \param[in, out] stp        pointer to cs_param_source_term_t structure
 * \param[in]      st_name    name of the source term
 * \param[in]      ml_id      id of the related to a cs_mesh_location_t struct.
 * \param[in]      type       type of source term
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      quad_type  type of quadrature rule to use
 * \param[in]      def_type   type of definition (by value, function...)
 * \param[in]      val        access to the definition of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_param_source_term_add(cs_param_source_term_t       *stp,
                         const char                   *st_name,
                         int                           ml_id,
                         cs_param_source_term_type_t   type,
                         cs_param_var_type_t           var_type,
                         cs_quadra_type_t              quad_type,
                         cs_param_def_type_t           def_type,
                         const void                   *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a reaction term
 *
 * \param[in] r_info     cs_param_reaction_t structure
 *
 * \return the name of the reaction term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_reaction_get_name(const cs_param_reaction_t   r_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of reaction term
 *
 * \param[in] r_info     set of parameters related to a reaction term
 *
 * \return the name associated with this type of reaction term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_reaction_get_type_name(cs_param_reaction_t  r_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a source term
 *
 * \param[in] st_info     cs_param_source_term_t structure
 *
 * \return the name of the source term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_source_term_get_name(const cs_param_source_term_t   st_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of type of a given source term structure
 *
 * \param[in] st_info     cs_param_source_term_t structure
 *
 * \return  the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_source_term_get_type_name(const cs_param_source_term_t   st_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of variable
 *
 * \param[in] type     cs_param_var_type_t
 *
 * \return the name associated to this type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_var_type_name(const cs_param_var_type_t   type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of definition
 *
 * \param[in] type     cs_param_def_type_t
 *
 * \return the name associated to this type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_def_type_name(const cs_param_def_type_t   type);

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
