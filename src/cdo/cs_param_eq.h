#ifndef __CS_PARAM_EQ_H__
#define __CS_PARAM_EQ_H__

/*============================================================================
 * Routines to handle the settings of a convection/diffusion equation
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo.h"
#include "cs_quadrature.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Boundary condition flag */
#define  CS_PARAM_EQ_UNSTEADY   (1 <<  0)  /*  1: unsteady term */
#define  CS_PARAM_EQ_CONVECTION (1 <<  1)  /*  2: convection term */
#define  CS_PARAM_EQ_DIFFUSION  (1 <<  2)  /*  4: diffusion term */
#define  CS_PARAM_EQ_SOURCETERM (1 <<  3)  /*  8: source term */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of equations managed by the solver */
typedef enum {

  CS_PARAM_EQ_TYPE_NONE,
  CS_PARAM_EQ_TYPE_SCAL,
  CS_PARAM_EQ_TYPE_VECT,
  CS_PARAM_EQ_TYPE_TENS,
  CS_PARAM_EQ_N_TYPES

} cs_param_eq_type_t;

/* Type of algorithm to get the solution of an equation */
typedef enum {

  CS_PARAM_EQ_ALGO_NONE,
  CS_PARAM_EQ_ALGO_ITSOL,   // Directly used an iterative solver
  CS_PARAM_EQ_ALGO_UZAWA,   // To solve sadle-point system
  CS_PARAM_EQ_ALGO_NEWTON,  // To solve non-linear system
  CS_PARAM_EQ_ALGO_PICARD,  // To solve non-linear system
  CS_PARAM_EQ_N_ALGOS

} cs_param_eq_algo_type_t;

/* Description of the algorithm used to solve an equation */
typedef struct {

  cs_param_eq_algo_type_t   type;

  int     n_iters;
  int     n_max_iters;
  int     n_cumulated_iters;
  int     n_max_cumulated_iters;

  double  eps;                    /* stopping criterion on accuracy */

} cs_param_eq_algo_t;

/* EQUATION STRUCTURE */
/* ================== */

/* Convection/Diffusion/Reaction + Term Source equation */
typedef struct {

  char *restrict       name;    /* Short description */

  cs_param_eq_type_t   type;    /* scalar, vector, tensor... */
  int                  iwarni;  /* Level of detail to output */

  /* Numerical settings */
  cs_space_scheme_t    space_scheme;

  /* Unsteady-Diffusion-Convection-Source term activated or not */
  int                  flag;

  /* Variable to solve (stored in a cs_field_t structure) */
  char *restrict       varname;
  int                  field_id;

  /* Unsteady term */
  cs_param_hodge_t     unsteady_hodge;
  bool                 is_multiplied_by_rho;  /* true or false */

  /* Diffusion parameters */
  cs_param_hodge_t     diffusion_hodge;

  /* Convection term (TODO) */

  /* Source term(s) */
  int                      n_source_terms;
  cs_param_source_term_t  *source_terms;

  /* Boundary conditions */
  cs_param_bc_t           *bc;

  /* High-level structure to manage/monitor the resolution of this equation */
  cs_param_eq_algo_t       algo_info;
  cs_param_itsol_t         itsol_info;

} cs_param_eq_t;


/* NAVSTO STRUCTURE */
/* ================ */

typedef enum {

  CS_PARAM_NAVSTO_FORM_NONE,
  CS_PARAM_NAVSTO_FORM_CURL,
  CS_PARAM_NAVSTO_FORM_CLASSIC,
  CS_PARAM_NAVSTO_N_FORMS

} cs_param_navsto_formulation_t;

typedef enum {

  CS_PARAM_NAVSTO_ALGO_NONE,
  CS_PARAM_NAVSTO_ALGO_COUPLED,
  CS_PARAM_NAVSTO_ALGO_SIMPLEC,
  CS_PARAM_NAVSTO_N_ALGOS

} cs_param_navsto_algo_t;

typedef struct {

  cs_param_navsto_formulation_t    formulation;
  cs_param_navsto_algo_t           algo;

  /* Steady/Reaction/Diffusion/Convection
     If diffusion is off  --> Euler equations
     If convection is off --> Stokes equations
   */
  cs_flag_t   momentum_flag;

  /* Fields (activated according to the formulation) */
  int    velocity_id;
  int    pressure_id;
  int    vorticity_id;

  /* Material property: rho and laminar viscosity */
  int    rho_id;
  int    lvisc_id;

  /* Source terms */
  int                      n_momemtum_source_terms;
  cs_param_source_term_t  *momentum_source_terms;

  int                      n_mass_source_terms;
  cs_param_source_term_t  *mass_source_terms;

  /* Boundary conditions */
  int                         n_bc_defs;
  cs_param_bc_def_t          *bc_defs;
  cs_param_bc_navsto_type_t   default_bc;

} cs_param_navsto_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the id related to an equation definition from its name
 *
 * \param[in]  ref_name    name of the property to find
 *
 * \return -1 if not found otherwise the associated id
 */
/*----------------------------------------------------------------------------*/

int
cs_param_eq_get_id_by_name(const char  *ref_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a cs_param_eq_t structure
 *
 * \param[in]  eq_id    id of the selected equation
 *
 * \return a pointer to the selected cs_param_eq_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_param_eq_t *
cs_param_eq_get_by_id(int   eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in] name             name of the material property
 * \param[in] varname          name of the variable associated to this equation
 * \param[in] type             type of equation (scalar, vector, tensor...)
 * \param[in] is_steady        add an unsteady term or not
 * \param[in] do_convection    add a convection term
 * \param[in] do_diffusion     add a diffusion term
 * \param[in] default_bc_type  type of boundary condition set by default
 *
 * \return  id associated to this equation
 */
/*----------------------------------------------------------------------------*/

int
cs_param_eq_add(const char               *name,
                const char               *varname,
                cs_param_eq_type_t        type,
                bool                      is_steady,
                bool                      do_convection,
                bool                      do_diffusion,
                cs_param_bc_basic_type_t  default_bc_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a material property to the diffusion term of an equation
 *         By default, a material property equal to the unity is set.
 *
 * \param[in]   eq_name    name of the equation to deal with
 * \param[in]   pty_name   name of the material property to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_set_diffusion_pty(const char   *eq_name,
                              const char   *pty_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Modify the level of warning (user function)
 *
 * \param[in]   name      name of the equation to deal with
 * \param[in]   iwarni    level of warning
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_set_warning_level(const char         *name,
                              int                 iwarni);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name associated to an equation from its id
 *
 * \param[in]   eq_id     id associated to a cs_param_eq_t structure
 *
 * \return the name of this equation
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_eq_get_name(int           eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Modify the scheme used to discretized in space an equation
 *
 * \param[in]   name      name of the equation to deal with
 * \param[in]   scheme    type of space scheme to use
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_set_space_scheme(const char         *name,
                             cs_space_scheme_t   scheme);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the numerical scheme used to discretize in space this
 *         equation from its id
 *
 * \param[in]   eq_id     id associated to a cs_param_eq_t structure
 *
 * \return the current space scheme used to discretize this equation
 */
/*----------------------------------------------------------------------------*/

cs_space_scheme_t
cs_param_eq_get_space_scheme(int    eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of equation associated to this equation from
 *         its id
 *
 * \param[in]   eq_id     id associated to a cs_param_eq_t structure
 *
 * \return the type of equation
 */
/*----------------------------------------------------------------------------*/

cs_param_eq_type_t
cs_param_eq_get_type(int    eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Modify the algorithm used to build the discrete Hodge operator
 *         related to the diffusion term
 *
 * \param[in]   name      name of the equation to deal with
 * \param[in]   scheme    type of space scheme to use
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_hodge_diffusion_set_algo(const char             *name,
                                     cs_param_hodge_algo_t   algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Modify the coefficient related to an algorithm used to build the
 *         discrete Hodge operator associated to the diffusion term
 *
 * \param[in]   name     name of the equation to deal with
 * \param[in]   coef     value of the coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_hodge_diffusion_set_coef(const char   *name,
                                     double        coef);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new boundary condition for a scalar
 *
 * \param[in]  eq_name        name of the equation
 * \param[in]  location_name  name of the location
 * \param[in]  bc_type        type of boundary condition
 * \param[in]  bc_val         value of the boundary condition
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_scalbc_by_val(const char                *eq_name,
                              const char                *location_name,
                              cs_param_bc_basic_type_t   bc_type,
                              double                     bc_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new boundary condition for a scalar using an analytic
 *         function
 *
 * \param[in]  eq_name        name of the equation
 * \param[in]  location_name  name of the location
 * \param[in]  bc_type        type of boundary condition
 * \param[in]  analytic       pointer to an analytic function
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_scalbc_by_analytic(const char                *eq_name,
                                   const char                *location_name,
                                   cs_param_bc_basic_type_t   bc_type,
                                   cs_analytic_func_t        *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a source term by value(s). This source term is added to the
 *         list of source terms associated to an equation
 *
 * \param[in]  eq_name   name of the equation
 * \param[in]  st_name   name of the source term (for log/post-processing)
 * \param[in]  ml_name   name of the mesh location
 * \param[in]  type      type of source term
 * \param[in]  get_imp   value(s) of the implicit part
 * \param[in]  get_exp   value(s) of the explicit part
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_source_term_by_val(const char                    *eq_name,
                                   const char                    *st_name,
                                   const char                    *ml_name,
                                   cs_param_source_term_type_t    type,
                                   cs_get_t                       get_imp,
                                   cs_get_t                       get_exp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a source term by a user-defined function. This source term
 *         is added to the list of source terms associated to an equation
 *
 * \param[in]  eq_name    name of the equation
 * \param[in]  st_name    name of the source term (for log/post-processing)
 * \param[in]  ml_name    name of the mesh location
 * \param[in]  type       type of source term
 * \param[in]  quad_type  quadrature rule
 * \param[in]  imp_func   pointer to a function related to the implicit part
 * \param[in]  exp_func   pointer to a function related to the explicit part
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_source_term_by_user(const char                   *eq_name,
                                    const char                   *st_name,
                                    const char                   *ml_name,
                                    cs_param_source_term_type_t   type,
                                    cs_quadra_type_t              quad_type,
                                    cs_user_func_t               *imp_func,
                                    cs_user_func_t               *exp_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a source term by an analytic function. This source term
 *         is added to the list of source terms associated to an equation
 *
 * \param[in]  eq_name    name of the equation
 * \param[in]  st_name    name of the source term (for log/post-processing)
 * \param[in]  ml_name    name of the mesh location
 * \param[in]  do_post    true or false
 * \param[in]  type       type of source term
 * \param[in]  quad_type  quadrature rule
 * \param[in]  imp_func   pointer to the function related to the implicit part
 * \param[in]  exp_func   pointer to the function related to the explicit part
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_source_term_by_analytic(const char                  *eq_name,
                                        const char                  *st_name,
                                        const char                  *ml_name,
                                        cs_param_source_term_type_t  type,
                                        cs_quadra_type_t             quad_type,
                                        cs_analytic_func_t          *imp_func,
                                        cs_analytic_func_t          *exp_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field related to a variable solved in an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all definitions of equations initialized during the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_free_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve high-level information on the setting
 *
 * \param[inout]   do_navsto     true or false
 * \param[inout]   n_cdo_eqs     number of additional equations using the CDO
 *                               kernel
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_get_info(bool        *do_navsto,
                     int         *n_cdo_eqs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Resume parameters of all conv./diff./source terms equations
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_resume_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_EQ_H__ */
