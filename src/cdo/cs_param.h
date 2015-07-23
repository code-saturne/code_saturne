#ifndef __CS_PARAM_H__
#define __CS_PARAM_H__

/*============================================================================
 * Manage the definition/setting of a computation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_cdo.h"
#include "cs_quadrature.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Tag to build parameter flag */
#define CS_PARAM_FLAG_UNIFORM  (1 <<  0)  //    1: uniform (in space)
#define CS_PARAM_FLAG_VERTEX   (1 <<  1)  //    2: on vertices
#define CS_PARAM_FLAG_EDGE     (1 <<  2)  //    4: on edges
#define CS_PARAM_FLAG_FACE     (1 <<  3)  //    8: on faces
#define CS_PARAM_FLAG_CELL     (1 <<  4)  //   16: on cells
#define CS_PARAM_FLAG_PRIMAL   (1 <<  5)  //   32: on primal mesh
#define CS_PARAM_FLAG_DUAL     (1 <<  6)  //   64: on dual mesh
#define CS_PARAM_FLAG_BORDER   (1 <<  7)  //  128: scalar-valued
#define CS_PARAM_FLAG_SCAL     (1 <<  8)  //  256: scalar-valued
#define CS_PARAM_FLAG_VECT     (1 <<  9)  //  512: vector-valued
#define CS_PARAM_FLAG_TENS     (1 << 10)  // 1024: tensor-valued
#define CS_PARAM_FLAG_SYMMET   (1 << 11)  // 2048: symmetric
#define CS_PARAM_FLAG_IMPLICIT (1 << 12)  // 4096: implicit
#define CS_PARAM_FLAG_EXPLICIT (1 << 13)  // 8192: explicit
#define CS_PARAM_FLAG_UNSTEADY (1 << 14)  //16384: unsteady

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef union {

  cs_get_t              get;       // definition by value
  cs_analytic_func_t   *analytic;  // definition by an analytic function
  cs_user_func_t       *user_func; // definition by an user-defined function

} cs_def_t;

typedef enum {

  CS_PARAM_DEF_NONE,
  CS_PARAM_DEF_BY_VALUE,
  CS_PARAM_DEF_BY_FIELD,
  CS_PARAM_DEF_BY_EVALUATOR,
  CS_PARAM_DEF_BY_ANALYTIC_FUNCTION,
  CS_PARAM_DEF_BY_USER_FUNCTION,
  CS_PARAM_DEF_BY_LAW,
  CS_PARAM_DEF_BY_FILE,
  CS_PARAM_N_DEF_TYPES

} cs_param_def_type_t;

/* Dimension of the variable to deal with */
typedef enum {

  CS_PARAM_VAR_NONE,
  CS_PARAM_VAR_SCAL,    // scalar variable (dim = 1)
  CS_PARAM_VAR_VECT,    // vector variable (dim = 3)
  CS_PARAM_VAR_SYMTENS, // symmetric tensor variable (dim = 6)
  CS_PARAM_VAR_TENS,    // tensor variable (dim = 9)
  CS_PARAM_N_VAR_TYPES

} cs_param_var_type_t;

/* MATERIAL PROPERTIES */
/* =================== */

/* Material property: parametrization for conductivity, diffusion coef.
   or viscosity */
typedef enum {

  CS_PARAM_PTY_NONE,
  CS_PARAM_PTY_ISO,
  CS_PARAM_PTY_ORTHO,
  CS_PARAM_PTY_ANISO,
  CS_PARAM_PTY_N_TYPES

} cs_param_pty_type_t;

/* Material parameter definitions */
typedef struct {

  char        *restrict name;
  cs_flag_t             flag;  /* Short descriptor (mask of bits) */

  int        post_freq; /* -1 (no post-processing), 0 (at the beginning)
                           otherwise every post_freq iteration(s) */
  int        field_id;  /* If post-processing is required (-1 if not used) */

  /* Material properties:
     Uniform material. The number of values to set depends on the material
     definition.
       - isotropic   = 1
       - orthotropic = 3
       - anisotropic = 9
     Not uniform or not steady material. Values have to be defined by a
     user function.
  */

  cs_param_pty_type_t  type;
  cs_param_def_type_t  def_type;
  cs_def_t             def;

} cs_param_pty_t;

/* DISCRETE HODGE OPERATORS */
/* ======================== */

typedef enum {

  CS_PARAM_HODGE_TYPE_NONE,
  CS_PARAM_HODGE_TYPE_VPCD,
  CS_PARAM_HODGE_TYPE_EPFD,
  CS_PARAM_HODGE_TYPE_FPED,
  CS_PARAM_HODGE_TYPE_EDFP,
  CS_PARAM_HODGE_TYPE_CPVD,
  CS_PARAM_HODGE_N_TYPES

} cs_param_hodge_type_t;

typedef enum {

  CS_PARAM_HODGE_ALGO_NONE,
  CS_PARAM_HODGE_ALGO_VORONOI,
  CS_PARAM_HODGE_ALGO_WBS,     // WBS: Whitney Barycentric Subdivision
  CS_PARAM_HODGE_ALGO_COST,    // COST: COnsistency & STabilization splitting
  CS_PARAM_HODGE_N_ALGOS

} cs_param_hodge_algo_t;

typedef struct {

  int                     pty_id;  /* id of the related material property */
  bool                    inv_pty; /* Definition based on material property
                                      or its inverse */

  cs_param_hodge_type_t   type;
  cs_param_hodge_algo_t   algo;
  double                  coef;    /* Value of the stabilization parameter
                                      if the COST algo. is used, other 0. */

} cs_param_hodge_t;

/* BOUNDARY CONDITIONS */
/* =================== */

/* Basic boundary conditions */
typedef enum {

  CS_PARAM_BC_BASIC_NONE,

  CS_PARAM_BC_HMG_DIRICHLET,
  CS_PARAM_BC_DIRICHLET,
  CS_PARAM_BC_HMG_NEUMANN,
  CS_PARAM_BC_NEUMANN,
  CS_PARAM_BC_ROBIN,

  CS_PARAM_BC_N_BASIC_TYPES

} cs_param_bc_basic_type_t;

/* Physic-driven boundary conditions */
typedef enum {

  CS_PARAM_BC_NAVSTO_NONE,

  CS_PARAM_BC_WALL,
  CS_PARAM_BC_MOVING_WALL,
  CS_PARAM_BC_INLET,
  CS_PARAM_BC_OUTLET,
  CS_PARAM_BC_SYMMETRY,

  /* Specific condition in CDO schemes for the Stokes problem in CURL
     formulation */
  CS_PARAM_BC_UNUT,   /* Normal and tangential component of velocity */
  CS_PARAM_BC_UNWT,   /* Normal compenent of velocity and viscosity times
                         tangential component of vorticity */
  CS_PARAM_BC_PRUT,   /* Tangential component of the velocity and pressure */

  CS_PARAM_BC_N_NAVSTO_TYPES

} cs_param_bc_navsto_type_t;

/* Values associated to the different ways to retrieve data */
typedef union {

  cs_param_bc_navsto_type_t   navsto;
  cs_param_bc_basic_type_t    basic;

} cs_param_bc_type_t;


/* Definition of the value of a boundary condition (BC)
   One or two values may be needed according to the type of BC.
   That's why two coefficients are used in the definition.
   For Dirichlet or Neumann BC for instance, only one coefficient is used.
*/

typedef struct {

  int                  loc_id;   // Id related to the list of border faces

  cs_param_bc_type_t   bc_type;
  cs_param_def_type_t  def_type; // Type of definition for a and b

  /* Access to the value related to the first coefficient and possibly
     the second one */
  cs_def_t             def_coef1;
  cs_def_t             def_coef2;

} cs_param_bc_def_t;

typedef struct {

  cs_param_bc_type_t   default_bc;
  bool                 strong_enforcement;
  double               penalty_coef; // TODO: preliminary step
  cs_quadra_type_t     quad_type;    // barycentric, higher, highest...

  int                  n_defs;
  cs_param_bc_def_t   *defs;

} cs_param_bc_t;

/* SOURCE TERMS */
/* ============ */

/* Types of source terms */
typedef enum {

  CS_PARAM_SOURCE_TERM_NONE,     // not set
  CS_PARAM_SOURCE_TERM_BASIC,    // in the right hand side
  CS_PARAM_SOURCE_TERM_IMPLICIT, // in the matrix to invert
  CS_PARAM_SOURCE_TERM_IMEX,     // implicit/explicit: two contributions
  CS_PARAM_SOURCE_TERM_MASS,     // specific treatment
  CS_PARAM_SOURCE_TERM_HEADLOSS, // specific treatment
  CS_PARAM_N_SOURCE_TERM_TYPES

} cs_param_source_term_type_t;

typedef struct {

  char  *restrict name;  /* short description of the source term */

  int            location_id;  /* list of cells */

  /* Specification related to the way of computing the source term */
  cs_param_source_term_type_t  type;      /* mass, head loss... */
  cs_param_var_type_t          var_type;  /* scalar, vector... */
  cs_param_def_type_t          def_type;  /* by value, by function... */
  cs_quadra_type_t             quad_type; /* barycentric, higher, highest */

  /* Two potential values (implicit and explicit) */
  cs_def_t                     imp_def;
  cs_def_t                     exp_def;

} cs_param_source_term_t;

/* ITERATIVE SOLVERS */
/* ================= */

typedef enum {

  CS_PARAM_PRECOND_NONE,
  CS_PARAM_PRECOND_DIAG,    // Diagonal preconditionning (also called Jacobi)
  CS_PARAM_PRECOND_POLY1,   // Neumann polynomial preconditionning (Order 1)
  CS_PARAM_PRECOND_SSOR,    // Symmetric Successive OverRelaxations
  CS_PARAM_PRECOND_ILU0,    // Incomplete LU factorization
  CS_PARAM_PRECOND_ICC0,    // Incomplete Cholesky factorization
  CS_PARAM_PRECOND_AMG,     // Algebraic MultiGrid
  CS_PARAM_PRECOND_AS,      // Additive Schwarz method
  CS_PARAM_N_PRECOND_TYPES

} cs_param_precond_type_t;

/* Type of iterative solver to use to inverse the linear system */
typedef enum {

  CS_PARAM_ITSOL_NONE,
  CS_PARAM_ITSOL_CG,
  CS_PARAM_ITSOL_BICG,
  CS_PARAM_ITSOL_GMRES,
  CS_PARAM_N_ITSOL_TYPES

} cs_param_itsol_type_t;

/* Description of the algorithm used to solve an equation */
typedef struct {

  cs_param_precond_type_t  precond;
  cs_param_itsol_type_t    solver;

  int                      n_max_iter;
  double                   eps;         // stopping criterion on accuracy

  int                      output_freq; // frequencdy of output into listing
  bool                     resid_normalized;

} cs_param_itsol_t;

/*============================================================================
 * Global variables
 *============================================================================*/

extern cs_flag_t  cs_glob_param_discretization_flag;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all structures related to properties
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_free_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a cs_param_pty_t structure from its id
 *
 * \param[in]  pty_id   id related to a property
 *
 * \return a pointer to a cs_param_pty_t
 */
/*----------------------------------------------------------------------------*/

cs_param_pty_t *
cs_param_pty_get(int  pty_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the id related to a property definition from its name
 *
 * \param[in]  ref_name    name of the property to find
 *
 * \return -1 if not found otherwise the associated id
 */
/*----------------------------------------------------------------------------*/

int
cs_param_pty_get_id_by_name(const char  *ref_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add by default several material properties
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_set_default(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and intialize a material property
 *
 * \param[in]  name        name of the material property
 * \param[in]  type        type of behavior of this material property
 * \param[in]  post_freq   -1 (no post-processing), 0 (at the beginning)
 *                         otherwise every post_freq iteration(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_add(const char             *name,
                 cs_param_pty_type_t     type,
                 int                     post_freq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field related to a material property
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_add_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a material property by value
 *
 * \param[in]  name      name of the material property
 * \param[in]  matval    value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_set_by_val(const char   *name,
                        cs_get_t      matval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a material property by an analytical function
 *
 * \param[in]  name       name of the material property
 * \param[in]  analytic   pointer to a cs_analytic_func_t
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_set_by_analytic_func(const char          *name,
                                  cs_analytic_func_t  *analytic_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query to know if the material property is uniform
 *
 * \param[in]    pty_id    id related to the material property to deal with
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_pty_is_uniform(int       pty_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of a material property from its id
 *
 * \param[in]    pty_id    id related to the material property to deal with
 *
 * \return  the name of the related property
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_pty_get_name(int            pty_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the 3x3 matrix related to a general material property.
 *         This value is computed at location (x,y,z) and time t.
 *
 * \param[in]    pty_id    id related to the material property to deal with
 * \param[in]    t         time at which we evaluate the material property
 * \param[in]    xyz       location at which  we evaluate the material property
 * \param[in]    invers    true or false
 * \param[inout] matval    pointer to the 3x3 matrix to return
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_get_val(int            pty_id,
                     cs_real_t      t,
                     cs_real_3_t    xyz,
                     bool           invers,
                     cs_real_33_t  *matval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Resume all the cs_param_pty_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_resume_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free structures dedicated to the definition of material properties
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_param_bc_t structure
 *
 * \param[in]  default_bc     default boundary condition
 * \param[in]  is_penalized   true/false
 *
 * \return a pointer to the new structure (free with cs_param_eq_t)
 */
/*----------------------------------------------------------------------------*/

cs_param_bc_t *
cs_param_bc_create(cs_param_bc_type_t  default_bc,
                   bool                is_penalized);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_param_bc_def_t structure
 *
 * \param[inout] bc_def     pointer to cs_param_bc_def_t struct. to set
 * \param[in]    loc_id     id related to a cs_mesh_location_t
 * \param[in]    bc_type    generic type of admissible boundary conditions
 * \param[in]    def_type   by value, function...
 * \param[in]    def_coef1  access to the value of the first coef
 * \param[in]    def_coef2  access to the value of the second coef (optional)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_bc_def_set(cs_param_bc_def_t         *bc_def,
                    int                        loc_id,
                    cs_param_bc_basic_type_t   bc_type,
                    cs_param_def_type_t        def_type,
                    cs_def_t                   def_coef1,
                    cs_def_t                   def_coef2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a source term. This source term is added to the list of
 *         source terms associated to an equation
 *
 * \param[inout] st         pointer to the cs_param_source_term_t struc. to set
 * \param[in]    st_name    name of the source term (for log/post-processing)
 * \param[in]    ml_id      id related to a cs_mesh_location_t
 * \param[in]    type       type of source term
 * \param[in]    var_type   type of variables (scalar, vector, tensor...)
 * \param[in]    quad_type  type of quadrature rule to use
 * \param[in]    def_type   type of definition (by value, function...)
 * \param[in]    imp_def    access to the definition of the implicit part
 * \param[in]    exp-def    access to the definition of the explicit part
 */
/*----------------------------------------------------------------------------*/

void
cs_param_source_term_add(cs_param_source_term_t       *st,
                         const char                   *st_name,
                         int                           ml_id,
                         cs_param_source_term_type_t   type,
                         cs_param_var_type_t           var_type,
                         cs_quadra_type_t              quad_type,
                         cs_param_def_type_t           def_type,
                         cs_def_t                      imp_def,
                         cs_def_t                      exp_def);

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
 * \brief   Get the name related to a source term
 *
 * \param[in] st_info     cs_param_source_term_t structure
 *
 * \return the name of the source term
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
 * \brief   Get the name of the preconditionner
 *
 * \param[in] precond     type of preconditionner
 *
 * \return the associated preconditionner name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_precond_name(cs_param_precond_type_t  precond);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_H__ */
