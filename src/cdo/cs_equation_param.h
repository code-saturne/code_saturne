#ifndef __CS_EQUATION_PARAM_H__
#define __CS_EQUATION_PARAM_H__

/*============================================================================
 * Header to handle specific settings related to a cs_equation_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param.h"
#include "cs_property.h"
#include "cs_advection_field.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Term flag */
#define CS_EQUATION_LOCKED        (1 <<  0)  //  1: modification not allowed
#define CS_EQUATION_UNSTEADY      (1 <<  1)  //  2: unsteady term
#define CS_EQUATION_CONVECTION    (1 <<  2)  //  4: convection term
#define CS_EQUATION_DIFFUSION     (1 <<  3)  //  8: diffusion term
#define CS_EQUATION_REACTION      (1 <<  4)  // 16: reaction term

/* Extra operations flag */
#define CS_EQUATION_POST_PECLET      (1 << 0) //  1: Export Peclet number
#define CS_EQUATION_POST_UPWIND_COEF (1 << 1) //  2: Export upwinding coef.

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of equations managed by the solver */
typedef enum {

  CS_EQUATION_TYPE_USER,         // User-defined equation
  CS_EQUATION_TYPE_GROUNDWATER,  // Equation specific to groundwater flows
  CS_EQUATION_TYPE_PREDEFINED,   // General predefined equation
  CS_EQUATION_N_TYPES

} cs_equation_type_t;

/* Type of algorithm to get the solution of an equation */
typedef enum {

  CS_EQUATION_ALGO_CS_ITSOL,    /* Used an iterative solver
                                   defined by Code_Saturne */
  CS_EQUATION_ALGO_PETSC_ITSOL, /* Used an iterative solver
                                   defined by PETSc */
  CS_EQUATION_ALGO_UZAWA,       // To solve sadle-point system
  CS_EQUATION_ALGO_NEWTON,      // To solve non-linear system
  CS_EQUATION_ALGO_PICARD,      // To solve non-linear system
  CS_EQUATION_N_ALGOS

} cs_equation_algo_type_t;

/* Description of the algorithm used to solve an equation */
typedef struct {

  cs_equation_algo_type_t   type;

  int     n_iters;
  int     n_max_iters;
  int     n_cumulated_iters;
  int     n_max_cumulated_iters;

  double  eps;                    /* stopping criterion on accuracy */

} cs_equation_algo_t;

/* Set of parameters to handle an unsteady convection-diffusion-reaction
   equation with term sources */
typedef struct {

  cs_equation_type_t    type;           /* predefined, user... */
  cs_param_var_type_t   var_type;       /* scalar, vector, tensor... */
  int                   verbosity;      /* Level of detail to output */
  int                   sles_verbosity; /* Level of detail to output for SLES */

  /* Unsteady-Diffusion-Convection-Source term activated or not */
  cs_flag_t             flag;

  /* Post-treatment */
  cs_flag_t             process_flag;  /* Type of post-treatment to do */

  /* Numerical settings */
  cs_space_scheme_t     space_scheme;
  int                   space_poly_degree; /* Max. degree of the polynomial
                                              basis */

  /* Boundary conditions */
  cs_param_bc_t        *bc;

  /* High-level structure to manage/monitor the resolution of this equation */
  cs_equation_algo_t    algo_info;
  cs_param_itsol_t      itsol_info;

  /* Unsteady term parametrization */
  cs_param_time_t       time_info;
  cs_param_hodge_t      time_hodge;
  cs_property_t        *time_property;

  /* Diffusion term parametrization */
  cs_param_hodge_t      diffusion_hodge;
  cs_property_t        *diffusion_property;

  /* Advection term parametrization */
  cs_param_advection_t  advection_info;
  cs_adv_field_t       *advection_field;

  /* Reaction term parametrization
     Belong to the left-hand and/or right-hand side according to the time scheme
  */
  cs_param_hodge_t      reaction_hodge;
  int                   n_reaction_terms;
  cs_property_t       **reaction_properties;
  cs_param_reaction_t  *reaction_info;

  /* Parametrization of source terms (Belong to the right-hand side) */
  int                   n_source_terms;
  cs_source_term_t     *source_terms;

} cs_equation_param_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_equation_param_t
 *
 * \param[in] type             type of equation
 * \param[in] var_type         type of variable (scalar, vector, tensor...)
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return a pointer to a new allocated cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_create(cs_equation_type_t     type,
                         cs_param_var_type_t    var_type,
                         cs_param_bc_type_t     default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_param_t
 *
 * \param[in] eqp          pointer to a cs_equation_param_t
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_free(cs_equation_param_t     *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_equation_param_t structure
 *
 * \param[in]  eqname   name of the related equation
 * \param[in]  eq       pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_summary(const char                 *eqname,
                          const cs_equation_param_t  *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize SLES structure for the resolution of the linear system
 *        according to the settings related to this equation
 *
 * \param[in]   eqname       pointer to an cs_equation_t structure
 * \param[in]   eqp          pointer to a cs_equation_param_t struct.
 * \param[in]   field_id     id of the cs_field_t struct. for this equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_init_sles(const char                 *eqname,
                            const cs_equation_param_t  *eqp,
                            int                         field_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_PARAM_H__ */
