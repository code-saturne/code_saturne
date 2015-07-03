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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_base.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_cdo.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_eq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local variables
 *============================================================================*/

cs_param_navsto_t  *_navsto_param = NULL;

int  cs_n_cdo_param_eqs = 0;
cs_param_eq_t  *cs_cdo_param_eqs = NULL;

/* Default initialization */
static cs_param_eq_algo_t _algo_default = {
  CS_PARAM_EQ_ALGO_ITSOL, // type of iterative solver
  0,                      // n_iters
  50,                     // max. number of iterations
  0,                      // n_cumulated_iters
  10000,                  // max. number of cumulated iterations
  1e-6                    // stopping criterion
};

static cs_param_itsol_t _itsol_default = {
  CS_PARAM_PRECOND_SSOR,  // preconditionner
  CS_PARAM_ITSOL_CG,      // iterative solver
  2500,                   // max. number of iterations
  1e-12,                  // stopping criterion on the accuracy
  150,                    // output frequency
  false                   // normalization of the residual (true or false)
};

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the id related to an equation name
 *         Output an error if this name is not already defined.
 *
 * \param[in]    eq_name        name of the equation
 * \param[inout] p_eq_id        pointer on the id of the equation
 */
/*----------------------------------------------------------------------------*/

static void
_check_eq_name(const char  *eq_name,
               int         *p_eq_id)
{
  /* Retrieve id of the equation and the location from their name */
  int  eq_id = cs_param_eq_get_id_by_name(eq_name);

  if (eq_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation name %s.\n"
                " This equation is not already defined.\n"), eq_name);

  /* Return pointers */
  *p_eq_id = eq_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the id related to an equation name
 *         Output an error if this name is not already defined.
 *
 * \param[in]    ml_name     name of the location
 * \param[inout] p_loc_id    pointer on the id of the location
 */
/*----------------------------------------------------------------------------*/

static void
_check_ml_name(const char   *ml_name,
                int         *p_loc_id)
{
  int  loc_id = cs_mesh_location_get_id_by_name(ml_name);

  if (loc_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);

  /* Return pointers */
  *p_loc_id = loc_id;
}

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
cs_param_eq_get_id_by_name(const char  *ref_name)
{
  int  i;

  int  eq_id = -1;
  int  reflen = strlen(ref_name);

  for (i = 0; i < cs_n_cdo_param_eqs; i++) {

    cs_param_eq_t  *eq = cs_cdo_param_eqs + i;
    int  len = strlen(eq->name);

    if (reflen == len) {
      if (strcmp(ref_name, eq->name) == 0) {
        eq_id = i;
        break;
      }
    }

  } /* Loops on mesh locations */

  return eq_id;
}

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
cs_param_eq_get_by_id(int   eq_id)
{
  if (eq_id < 0 || eq_id >= cs_n_cdo_param_eqs)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation id (= %d). Stop execution.\n"), eq_id);

  return  cs_cdo_param_eqs + eq_id;
}

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
                _Bool                     is_steady,
                _Bool                     do_convection,
                _Bool                     do_diffusion,
                cs_param_bc_basic_type_t  default_bc_basic_type)
{
  int  len = strlen(name)+1;
  int  eq_id = cs_param_eq_get_id_by_name(name);
  cs_param_bc_type_t  default_bc_type;
  cs_param_eq_t  *eq = NULL;

  if (eq_id > -1) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_(" An existing equation has already the same name %s.\n"
                 " Stop adding the equation.\n"), name);
    return  eq_id;
  }

  eq_id = cs_n_cdo_param_eqs;
  cs_n_cdo_param_eqs++;
  BFT_REALLOC(cs_cdo_param_eqs, cs_n_cdo_param_eqs, cs_param_eq_t);
  eq = cs_cdo_param_eqs + eq_id;

  /* Initialize the equation structure by default */
  BFT_MALLOC(eq->name, len, char);
  strncpy(eq->name, name, len);

  eq->type = type;
  eq->iwarni = 0;
  eq->space_scheme = CS_SPACE_SCHEME_CDOVB;
  eq->field_id = -1;  // field is created when all user-defined data are set
  eq->is_multiplied_by_rho = true;
  eq->n_source_terms = 0;
  eq->source_terms = NULL;

  /* Build the equation flag */
  eq->flag = 0;
  if (!is_steady) {
    eq->flag |= CS_PARAM_EQ_UNSTEADY;

    eq->is_multiplied_by_rho = true;

    /* Default initialization in accordance with the default settings */
    eq->unsteady_hodge.pty_id = 0;      // Unity (default property)
    eq->unsteady_hodge.inv_pty = false; // inverse property ?
    eq->unsteady_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
    eq->unsteady_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
  }

  if (do_convection) {
    eq->flag |= CS_PARAM_EQ_CONVECTION;
  }

  if (do_diffusion) {
    eq->flag |= CS_PARAM_EQ_DIFFUSION;

    /* Default: DGA Hodge operator */
    eq->diffusion_hodge.pty_id = 0;      // Unity (default property)
    eq->diffusion_hodge.inv_pty = false; // inverse property ?
    eq->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    eq->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    eq->diffusion_hodge.coef = 1./3.;
  }

  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No variable name associated to equation id %d\n"
                " Check your initialization.\n"), eq_id);
  else {
    len = strlen(varname)+1;
    BFT_MALLOC(eq->varname, len, char);
    strncpy(eq->varname, varname, len);
  }

  default_bc_type.basic = default_bc_basic_type;
  eq->bc = cs_param_bc_create(default_bc_type,
                              false);           // penalized BCs ?

  eq->algo_info = _algo_default;
  eq->itsol_info = _itsol_default;

  return eq_id;
}

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
                              const char   *pty_name)
{
  int  eq_id = cs_param_eq_get_id_by_name(eq_name);
  int  pty_id = cs_param_pty_get_id_by_name(pty_name);

  /* Sanity checks */
  if (eq_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation name %s.\n"
                " This equation is not already defined.\n"), eq_name);

  if (pty_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid material property name %s.\n"
                " This material property is not already defined.\n"), pty_name);

  /* Associate the material property */
  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  eq->diffusion_hodge.pty_id = pty_id;
}

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
                              int                 iwarni)
{
  int  eq_id = cs_param_eq_get_id_by_name(name);

  /* Sanity checks */
  if (eq_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation name %s.\n"
                " This equation is not already defined.\n"), name);

  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  eq->iwarni = iwarni;
}

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
cs_param_eq_get_name(int           eq_id)
{
  if (eq_id < 0 || eq_id >= cs_n_cdo_param_eqs)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation id (= %d)\n"), eq_id);

  return cs_cdo_param_eqs[eq_id].name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Modify the scheme used to discretize in space an equation
 *         User function
 *
 * \param[in]   name      name of the equation to deal with
 * \param[in]   scheme    type of space scheme to use
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_set_space_scheme(const char         *name,
                             cs_space_scheme_t   scheme)
{
  int  eq_id = cs_param_eq_get_id_by_name(name);
  cs_param_eq_t  *eq = NULL;

  if (eq_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation name %s.\n"
                " This equation is not already defined.\n"), name);

  eq = cs_cdo_param_eqs + eq_id;
  eq->space_scheme = scheme;

  switch(scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    eq->unsteady_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
    eq->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    break;

  case CS_SPACE_SCHEME_CDOFB:
    eq->unsteady_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
    eq->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid space discretization scheme."));
    break;
  }

}

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
cs_param_eq_get_space_scheme(int                 eq_id)
{
  if (eq_id < 0 || eq_id >= cs_n_cdo_param_eqs)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation id (= %d)\n"), eq_id);

  return cs_cdo_param_eqs[eq_id].space_scheme;
}

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
cs_param_eq_get_type(int    eq_id)
{
  if (eq_id < 0 || eq_id >= cs_n_cdo_param_eqs)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation id (= %d)\n"), eq_id);

  return cs_cdo_param_eqs[eq_id].type;
}

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
                                     cs_param_hodge_algo_t   algo)
{
  int  eq_id = cs_param_eq_get_id_by_name(name);

  /* Sanity checks */
  if (eq_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation name %s.\n"
                " This equation is not already defined.\n"), name);

  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  eq->diffusion_hodge.algo = algo;
}

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
                                     double        coef)
{
  int  eq_id = cs_param_eq_get_id_by_name(name);

  /* Sanity checks */
  if (eq_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation name %s.\n"
                " This equation is not already defined.\n"), name);

  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  if (eq->diffusion_hodge.algo != CS_PARAM_HODGE_ALGO_COST) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_(" The algorithm used to build the discrete Hodge operator\n"
                 " related to the diffusion term doesn't need a coefficient.\n"
                 " Please check your settings.\n"));
  }

  eq->diffusion_hodge.coef = coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new boundary condition for a scalar equation
 *
 * \param[in]  eq_name    name of the equation
 * \param[in]  ml_name    name of the mesh location
 * \param[in]  bc_type    type of boundary condition
 * \param[in]  bc_val     value of the boundary condition
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_scalbc_by_val(const char               *eq_name,
                              const char               *ml_name,
                              cs_param_bc_basic_type_t  bc_type,
                              double                    bc_val)
{
  int  eq_id, ml_id, def_id;
  cs_def_t  def_coef1, def_coef2;

  /* Retrieve related ids */
  _check_eq_name(eq_name, &eq_id);
  _check_ml_name(ml_name, &ml_id);

  cs_param_eq_t *eq = cs_cdo_param_eqs + eq_id;
  cs_param_bc_t  *bc = eq->bc;

  /* Sanity checks */
  assert(bc != NULL);
  assert(eq->type == CS_PARAM_EQ_TYPE_SCAL);

  /* Add a new definition */
  def_id = bc->n_defs;
  bc->n_defs += 1;
  BFT_REALLOC(bc->defs, bc->n_defs, cs_param_bc_def_t);

  def_coef1.get.val = bc_val;
  def_coef2.get.val = 0.; // Not used (default value).

  cs_param_bc_def_set(bc->defs + def_id,
                      ml_id,
                      bc_type,
                      CS_PARAM_DEF_BY_VALUE,
                      def_coef1, def_coef2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new boundary condition for a scalar using a user-defined
 *         function
 *
 * \param[in]  eq_name        name of the equation
 * \param[in]  location_name  name of the mesh location
 * \param[in]  bc_type        type of boundary condition
 * \param[in]  analytic       pointer to an analytic function
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_scalbc_by_analytic(const char                *eq_name,
                                   const char                *ml_name,
                                   cs_param_bc_basic_type_t   bc_type,
                                   cs_analytic_func_t        *analytic)
{
  int  eq_id, ml_id, def_id;
  cs_def_t  def_coef1, def_coef2;

  /* Retrieve related ids */
  _check_eq_name(eq_name, &eq_id);
  _check_ml_name(ml_name, &ml_id);

  cs_param_eq_t *eq = cs_cdo_param_eqs + eq_id;
  cs_param_bc_t  *bc = eq->bc;

  /* Sanity checks */
  assert(bc != NULL);
  assert(eq->type == CS_PARAM_EQ_TYPE_SCAL);

  /* Add a new definition */
  def_id = bc->n_defs;
  bc->n_defs += 1;
  BFT_REALLOC(bc->defs, bc->n_defs, cs_param_bc_def_t);

  def_coef1.analytic = analytic;
  def_coef2.get.val = 0; // Not used (default value).

  cs_param_bc_def_set(bc->defs + def_id,
                      ml_id,
                      bc_type,
                      CS_PARAM_DEF_BY_ANALYTIC_FUNCTION,
                      def_coef1, def_coef2);
}

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
                                   cs_get_t                       get_exp)
{
  int  k, l, eq_id, ml_id, st_id;
  cs_param_var_type_t  var_type;
  cs_def_t  imp_def, exp_def;

  /* Retrieve related id */
  _check_eq_name(eq_name, &eq_id);
  _check_ml_name(ml_name, &ml_id);

  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  st_id = eq->n_source_terms;
  eq->n_source_terms += 1;
  BFT_REALLOC(eq->source_terms, eq->n_source_terms, cs_param_source_term_t);

  switch (eq->type) {
  case CS_PARAM_EQ_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    imp_def.get.val = get_imp.val;
    exp_def.get.val = get_exp.val;
    break;

  case CS_PARAM_EQ_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    for (k = 0; k < 3; k++) {
      imp_def.get.vect[k] = get_imp.vect[k];
      exp_def.get.vect[k] = get_exp.vect[k];
    }
    break;

  case CS_PARAM_EQ_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    for (k = 0; k < 3; k++) {
      for (l = 0; l < 3; l++) {
        imp_def.get.tens[k][l] = get_imp.tens[k][l];
        exp_def.get.tens[k][l] = get_exp.tens[k][l];
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation. This case is not managed.\n"));
    break;
  }

  cs_param_source_term_add(eq->source_terms + st_id,
                           st_name,
                           ml_id,
                           type,
                           var_type,
                           CS_QUADRATURE_BARY,
                           CS_PARAM_DEF_BY_VALUE,
                           imp_def,
                           exp_def);
}

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
                                    cs_user_func_t               *exp_func)
{
  int  eq_id, ml_id, st_id;
  cs_param_var_type_t  var_type;
  cs_def_t  imp_def, exp_def;

  /* Retrieve related id */
  _check_eq_name(eq_name, &eq_id);
  _check_ml_name(ml_name, &ml_id);

  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  st_id = eq->n_source_terms;
  eq->n_source_terms += 1;
  BFT_REALLOC(eq->source_terms, eq->n_source_terms, cs_param_source_term_t);

  switch (eq->type) {
  case CS_PARAM_EQ_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    break;
  case CS_PARAM_EQ_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    break;
  case CS_PARAM_EQ_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation. This case is not managed.\n"));
    break;
  }

  imp_def.user_func = imp_func;
  exp_def.user_func = exp_func;

  cs_param_source_term_add(eq->source_terms + st_id,
                           st_name,
                           ml_id,
                           type,
                           var_type,
                           quad_type,
                           CS_PARAM_DEF_BY_USER_FUNCTION,
                           imp_def,
                           exp_def);
}

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
                                        cs_analytic_func_t          *exp_func)
{
  int  eq_id, ml_id, st_id;
  cs_param_var_type_t  var_type;
  cs_def_t  imp_def, exp_def;

  /* Retrieve related id */
  _check_eq_name(eq_name, &eq_id);
  _check_ml_name(ml_name, &ml_id);

  cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

  st_id = eq->n_source_terms;
  eq->n_source_terms += 1;
  BFT_REALLOC(eq->source_terms, eq->n_source_terms, cs_param_source_term_t);

  switch (eq->type) {
  case CS_PARAM_EQ_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    break;
  case CS_PARAM_EQ_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    break;
  case CS_PARAM_EQ_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation. This case is not managed.\n"));
    break;
  }

  imp_def.analytic = imp_func;
  exp_def.analytic = exp_func;

  cs_param_source_term_add(eq->source_terms + st_id,
                           st_name,
                           ml_id,
                           type,
                           var_type,
                           quad_type,
                           CS_PARAM_DEF_BY_ANALYTIC_FUNCTION,
                           imp_def,
                           exp_def);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field related to a variable solved in an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_add_fields(void)
{
  int  eq_id, dim, location_id;

  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

  for (eq_id = 0; eq_id < cs_n_cdo_param_eqs; eq_id++) {

    cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;
    _Bool has_previous = (eq->flag & CS_PARAM_EQ_UNSTEADY) ? true : false;

    /* Define dim */
    switch (eq->type) {
    case CS_PARAM_EQ_TYPE_SCAL:
      dim = 1;
      break;
    case CS_PARAM_EQ_TYPE_VECT:
      dim = 3;
      break;
    case CS_PARAM_EQ_TYPE_TENS:
      dim = 9;
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Type of equation for eq. %s is incompatible with the"
                  " creation of field.\n"
                  " Stop adding field from CDO user equations.\n"), eq->name);
    }

    /* Define mesh_location_id */
    switch (eq->space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
      location_id = cs_mesh_location_get_id_by_name(N_("vertices"));
      break;
    case CS_SPACE_SCHEME_CDOFB:
      location_id = cs_mesh_location_get_id_by_name(N_("cells"));
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Space scheme for eq. %s is incompatible with the"
                  " creation of field.\n"
                  " Stop adding field from CDO user equations.\n"), eq->name);
    }

    if (location_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid mesh location id (= -1) for the current field\n"));

    cs_field_t  *fld = cs_field_create(eq->varname,
                                       field_mask,
                                       location_id,
                                       dim,
                                       true,          // interleave
                                       has_previous);

    /* Allocate and initializa values */
    cs_field_allocate_values(fld);

    eq->field_id = cs_field_id_by_name(eq->varname);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all definitions of equations initialized during the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_free_all(void)
{
  int  i, eq_id;

  /* Free navsto structure (TODO) */

  /* Other equations */
  for (eq_id = 0; eq_id < cs_n_cdo_param_eqs; eq_id++) {

    cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;

    BFT_FREE(eq->name);
    BFT_FREE(eq->varname);

    if (eq->bc != NULL) { // Boundary conditions
      if (eq->bc->n_defs > 0)
        BFT_FREE(eq->bc->defs);
      BFT_FREE(eq->bc);
      eq->bc = NULL;
    }

    if (eq->n_source_terms > 0) { // Source terms
      for (i = 0; i< eq->n_source_terms; i++)
        BFT_FREE(eq->source_terms[i].name);
      BFT_FREE(eq->source_terms);
    }

  } /* Loop on equations */

  BFT_FREE(cs_cdo_param_eqs);
  cs_cdo_param_eqs = NULL;
  cs_n_cdo_param_eqs = 0;
}

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
cs_param_eq_get_info(_Bool       *do_navsto,
                     int         *n_cdo_eqs)
{
  *do_navsto = false;
  if (_navsto_param != NULL)
    *do_navsto = true;

  *n_cdo_eqs = cs_n_cdo_param_eqs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Resume parameters of all conv./diff./source terms equations
 */
/*----------------------------------------------------------------------------*/

void
cs_param_eq_resume_all(void)
{
  int  eq_id;

  if (cs_n_cdo_param_eqs == 0)
    bft_printf("  -msg- No additional equation is defined.\n");

  for (eq_id = 0; eq_id < cs_n_cdo_param_eqs; eq_id++) {

    const cs_param_eq_t  *eq = cs_cdo_param_eqs + eq_id;
    const cs_param_itsol_t   itsol = eq->itsol_info;

    bft_printf("\n");
    bft_printf(lsepline);
    bft_printf("  Resume settings for %s eq. (variable %s)\n",
               eq->name, eq->varname);
    bft_printf(lsepline);

    _Bool  unsteady = (eq->flag & CS_PARAM_EQ_UNSTEADY) ? true : false;
    _Bool  convection = (eq->flag & CS_PARAM_EQ_CONVECTION) ? true : false;
    _Bool  diffusion = (eq->flag & CS_PARAM_EQ_DIFFUSION) ? true : false;
    _Bool  source_term = (eq->n_source_terms > 0) ? true : false;

    bft_printf("  <Equation>  unsteady [%s], convection [%s], diffusion [%s],"
               " source term [%s]\n",
               cs_base_strtf(unsteady), cs_base_strtf(convection),
               cs_base_strtf(diffusion), cs_base_strtf(source_term));

    /* if (eqp->convection) { */
    /*   bft_printf("    -> Given convection field [%s]", */
    /*              cs_base_strtf(eqp->given_field)); */
    /*   if (eqp->given_field) */
    /*     bft_printf(" field: %p\n", (const void *)eqp->convfield); */
    /* } */

    if (eq->flag & CS_PARAM_EQ_DIFFUSION) {

      const cs_param_hodge_t  h_info = eq->diffusion_hodge;

      bft_printf("\t<Equation/Diffusion term>\n");
      bft_printf("\t\t--> Property related to the diffusion term: %s\n",
                 cs_param_pty_get_name(h_info.pty_id));

      if (eq->iwarni > 0) {
        bft_printf("\t\t--> Hodge operator: %s / %s\n",
                   cs_param_hodge_get_type_name(h_info),
                   cs_param_hodge_get_algo_name(h_info));
        bft_printf("\t\t--> Inversion of the material property: %s\n",
                   cs_base_strtf(h_info.inv_pty));
        if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
          bft_printf("\t\t--> Coefficient value for COST algo: %.3e\n",
                     h_info.coef);
      }

    } // Diffusion term

    if (eq->n_source_terms > 0) {

      int  s_id;

      bft_printf("\t<Equation/Source terms>\n");
      for (s_id = 0; s_id < eq->n_source_terms; s_id++) {

        cs_param_source_term_t  st_info = eq->source_terms[s_id];

        bft_printf("\t\t--> Source term label: %s\n",
                   cs_param_source_term_get_name(st_info));
        bft_printf("\t\t--> Related mesh location: %s\n",
                   cs_mesh_location_get_name(st_info.location_id));
        bft_printf("\t\t--> Type: %s; Variable type: %s; Definition type: %s\n",
                   cs_param_source_term_get_type_name(st_info),
                   cs_param_get_var_type_name(st_info.var_type),
                   cs_param_get_def_type_name(st_info.def_type));
        if (eq->iwarni > 0)
          bft_printf("\t\t--> Quadrature type: %s\n",
                     cs_quadrature_get_type_name(st_info.quad_type));

      } // Loop on source terms

    } // Source terms

    /* /\* Boundary definitions *\/ */
    /* bft_printf("\n"); */
    /* if (eqp->penalized_bc) */
    /*   bft_printf("  Penalization of boundary condition [True]; coef = %5.3e\n", */
    /*              eqp->penality_coef); */
    /* else */
    /*   bft_printf("  Penalization of boundary condition [False];\n"); */

    /* cs_param_bc_resume(eqp->n_bc_defs, eqp->bc_defs, eqp->default_bc); */

    /* Iterative solver information */
    bft_printf("\n  <Iterative Solver Parameters>\n");
    bft_printf(" -sla- Solver.MaxIter     %d\n", itsol.n_max_iter);
    bft_printf(" -sla- Solver.Name        %s\n",
               cs_param_get_solver_name(itsol.solver));
    bft_printf(" -sla- Solver.Precond     %s\n",
               cs_param_get_precond_name(itsol.precond));
    bft_printf(" -sla- Solver.Eps        % -10.6e\n", itsol.eps);
    bft_printf(" -sla- Solver.Normalized  %s\n",
               cs_base_strtf(itsol.resid_normalized));

    /* if (eqp->saddle_solver.type != CS_SADDLE_ALGO_NONE) */
    /*   cs_saddle_info_resume(eqp->saddle_solver); */

  } /* Loop on equations */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
