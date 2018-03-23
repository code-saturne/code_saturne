/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_hho_scaleq.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_range_set.h"
#include "cs_sles.h"
#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_EQUATION_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

static int  _n_equations = 0;
static int  _n_predef_equations = 0;
static int  _n_user_equations = 0;
static cs_equation_t  **_equations = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_eq[] =
  N_(" Stop setting an empty cs_equation_t structure.\n"
     " Please check your settings.\n");

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values for the variable related to an equation
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_initialize_field_from_ic(cs_equation_t  *eq)
{
  assert(eq != NULL);
  const cs_equation_param_t  *eqp = eq->param;

  /* Retrieve the associated field */
  cs_field_t  *field = cs_field_by_id(eq->field_id);
  cs_real_t  *values = field->val;

  cs_flag_t  dof_flag = 0;
  switch (eqp->dim) {
  case 1: /* Scalar-valued */
    dof_flag |= CS_FLAG_SCALAR;
    break;
  case 3: /* Vector-valued */
    dof_flag |= CS_FLAG_VECTOR;
    break;
  case 9: /* Tensor-valued */
    dof_flag |= CS_FLAG_TENSOR;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Incompatible type of variable for equation %s."), eq->name);
     break;
  }

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB) {

    cs_flag_t  v_flag = dof_flag | cs_flag_primal_vtx;

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(v_flag, def, values);
        break;

      case CS_XDEF_BY_QOV:
        cs_evaluate_potential_by_qov(v_flag, def, values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_potential_by_analytic(v_flag, def, values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

  } // VB or VCB schemes --> initialize on vertices

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB ||
      eqp->space_scheme == CS_SPACE_SCHEME_HHO_P0) {

    cs_flag_t  f_flag = dof_flag | cs_flag_primal_face;
    cs_real_t  *f_values = eq->get_extra_values(eq->builder);
    assert(f_values != NULL);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];

      /* Initialize face-based array */
      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(f_flag, def, f_values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_potential_by_analytic(f_flag, def, f_values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

  } // FB schemes --> initialize on faces

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB ||
      eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB ||
      eqp->space_scheme == CS_SPACE_SCHEME_HHO_P0) {

    /* TODO: HHO */

    /* Initialize cell-based array */
    cs_flag_t  c_flag = dof_flag | cs_flag_primal_cell;
    cs_real_t  *c_values = values;
    if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
      c_values = eq->get_extra_values(eq->scheme_context);
    assert(c_values != NULL);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];

      /* Initialize cell-based array */
      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(c_flag, def, c_values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_potential_by_analytic(c_flag, def, c_values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

  } // FB or VCB schemes --> initialize on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_to_cast   pointer to generic builder structure
 * \param[in, out] p_x          pointer of pointer to the solution array
 * \param[in, out] p_rhs        pointer of pointer to the RHS array
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_vb_solving(void              *eq_to_cast,
                    cs_real_t         *p_x[],
                    cs_real_t         *p_rhs[])
{
  cs_equation_t  *eq = (cs_equation_t  *)eq_to_cast;
  const cs_field_t  *fld = cs_field_by_id(eq->field_id);
  const int  eq_dim = fld->dim;

  cs_real_t  *x = NULL, *b = NULL;

  BFT_MALLOC(x, CS_MAX(eq->n_sles_scatter_elts,
                       cs_matrix_get_n_columns(eq->matrix)), cs_real_t);

  /* x and b are a "gathered" view of field->val and eq->rhs respectively
     through the range set operation.
     Their size is equal to n_sles_gather_elts <= n_sles_scatter_elts */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    /* Compact numbering to fit the algebraic decomposition */
    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        fld->val,             // in: size = n_sles_scatter_elts
                        x);                   //out: size = n_sles_gather_elts

    /* The right-hand side stems from a cellwise building on this rank.
       Other contributions from distant ranks may contribute to an element
       owned by the local rank */
    BFT_MALLOC(b, eq->n_sles_scatter_elts, cs_real_t);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      b[i] = eq->rhs[i];
#else
    memcpy(b, eq->rhs, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    cs_interface_set_sum(eq->rset->ifs,
                         eq->n_sles_scatter_elts, eq_dim, false, CS_REAL_TYPE,
                         b);

    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        b,                    // in: size = n_sles_scatter_elts
                        b);                   //out: size = n_sles_gather_elts

  }
  else { /* Serial mode *** without periodicity *** */

    assert(eq->n_sles_gather_elts == eq->n_sles_scatter_elts);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      x[i] = fld->val[i];
#else
    memcpy(x, fld->val, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    /* Nothing to do for the right-hand side */
    b = eq->rhs;

  }

  /* Return pointers */
  *p_x = x;
  *p_rhs = b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_to_cast   pointer to generic builder structure
 * \param[in, out] p_x          pointer of pointer to the solution array
 * \param[in, out] p_rhs        pointer of pointer to the RHS array
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_fb_solving(void              *eq_to_cast,
                    cs_real_t         *p_x[],
                    cs_real_t         *p_rhs[])
{
  cs_equation_t  *eq = (cs_equation_t  *)eq_to_cast;
  const cs_field_t  *fld = cs_field_by_id(eq->field_id);
  const cs_real_t  *f_values = eq->get_extra_values(eq->scheme_context);
  const int  eq_dim = fld->dim;

  /* Sanity check */
  assert(f_values != NULL);

  cs_real_t  *x = NULL, *b = NULL;
  BFT_MALLOC(x, CS_MAX(eq->n_sles_scatter_elts,
                       cs_matrix_get_n_columns(eq->matrix)), cs_real_t);

  /* x and b are a "gathered" view of field->val and eq->rhs respectively
     through the range set operation.
     Their size is equal to n_sles_gather_elts <= n_sles_scatter_elts */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    /* Compact numbering to fit the algebraic decomposition */
    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        f_values,             // in: size = n_sles_scatter_elts
                        x);                   //out: size = n_sles_gather_elts

    /* The right-hand side stems from a cellwise building on this rank.
       Other contributions from distant ranks may contribute to an element
       owned by the local rank */
    BFT_MALLOC(b, eq->n_sles_scatter_elts, cs_real_t);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      b[i] = eq->rhs[i];
#else
    memcpy(b, eq->rhs, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    cs_interface_set_sum(eq->rset->ifs,
                         eq->n_sles_scatter_elts, eq_dim, false, CS_REAL_TYPE,
                         b);

    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        b,                    // in: size = n_sles_scatter_elts
                        b);                   //out: size = n_sles_gather_elts

  }
  else { /* Serial mode *** without periodicity *** */

    assert(eq->n_sles_gather_elts == eq->n_sles_scatter_elts);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      x[i] = f_values[i];
#else
    memcpy(x, f_values, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    /* Nothing to do for the right-hand side */
    b = eq->rhs;

  }

  /* Return pointers */
  *p_x = x;
  *p_rhs = b;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of equations
 *
 * \return the current number of cs_equation_t structure allocated
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_n_equations(void)
{
  return _n_equations;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure with name eqname
 *         Return NULL if not find
 *
 * \param[in]  eqname    name of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_name(const char    *eqname)
{
  cs_equation_t  *eq = NULL;
  if (eqname == NULL)
    return eq;

  size_t  len_in = strlen(eqname);
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *_eq = _equations[i];
    if (strlen(_eq->name) == len_in)
      if (strcmp(eqname, _eq->name) == 0)
        return _eq;

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure associated to a
 *         cs_equation_t structure thanks to the equation name
 *
 * \param[in]  eqname       name of the equation
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_by_name(const char    *eqname)
{
  if (eqname == NULL)
    return NULL;

  else {

    cs_equation_t  *eq = cs_equation_by_name(eqname);

    if (eq == NULL)
      return NULL;
    else
      return eq->param;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure associated to a
 *         cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_get_param(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure with name eqname
 *         Return NULL if not find
 *
 * \param[in]  eq_id    id of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_id(int   eq_id)
{
  if (eq_id < 0 || eq_id > _n_equations - 1)
    return NULL;

  return _equations[eq_id];
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the name related to the given cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a name or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const char *
cs_equation_get_name(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the id number related to the given cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return an id (0 ... n-1) or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_id(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return -1;
  else
    return eq->id;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the field structure associated to a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_equation_get_field(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return cs_field_by_id(eq->field_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the flag associated to an equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a flag (cs_flag_t type)
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_equation_get_flag(const cs_equation_t    *eq)
{
  cs_flag_t  ret_flag = 0;

  if (eq == NULL)
    return ret_flag;

  ret_flag = eq->param->flag;

  return ret_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_builder_t structure associated to a
 *         cs_equation_t structure. Only for an advanced usage.
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_equation_builder_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_builder_t *
cs_equation_get_builder(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->builder;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to a structure useful to handle low-level
 *         operations for the given equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a structure to cast on-the-fly or NULL if not found
 */
/*----------------------------------------------------------------------------*/

void *
cs_equation_get_scheme_context(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->scheme_context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         diffusion term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_diffusion_property(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else if (eq->param == NULL)
    return NULL;
  else
    return eq->param->diffusion_property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         unsteady term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_time_property(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else if (eq->param == NULL)
    return NULL;
  else
    return eq->param->time_property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         reaction term called r_name and related to this equation
 *
 *
 * \param[in]  eq            pointer to a cs_equation_t structure
 * \param[in]  reaction_id   id related to this reaction term
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_reaction_property(const cs_equation_t    *eq,
                                  const int               reaction_id)
{
  if (eq == NULL)
    return NULL;

  const cs_equation_param_t  *eqp = eq->param;
  if (reaction_id < 0 || reaction_id > eqp->n_reaction_terms - 1)
    return NULL;

  return eqp->reaction_properties[reaction_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return true is the given equation is steady otherwise false
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_is_steady(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return true;
  if (eq->param == NULL)
    return true;

  if (eq->param->flag & CS_EQUATION_UNSTEADY)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of numerical scheme used for the discretization in
 *         space
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_param_space_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_param_space_scheme_t
cs_equation_get_space_scheme(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_SPACE_N_SCHEMES;
  else if (eq->param == NULL)
    return CS_SPACE_N_SCHEMES;
  else
    return eq->param->space_scheme;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the max. degree used in the polynomial basis for the space
 *         discretization
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the polynomial order
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_space_poly_degree(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return -1;
  else if (eq->param == NULL)
    return -1;
  else
    return eq->param->space_poly_degree;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the dimension of the variable solved by this equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  an integer corresponding to the dimension of the variable
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_var_dim(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return 0;
  else if (eq->param == NULL)
    return 0;
  else
    return eq->param->dim;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of equation for the given equation structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the type of the given equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_type_t
cs_equation_get_type(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_EQUATION_N_TYPES;
  else if (eq->param == NULL)
    return CS_EQUATION_N_TYPES;
  else
    return eq->param->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation structure and set a first set of parameters
 *
 * \param[in] eqname        name of the equation
 * \param[in] varname       name of the variable associated to this equation
 * \param[in] eqtype        type of equation (user, predefined...)
 * \param[in] dim           dimension of the unknow attached to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add(const char            *eqname,
                const char            *varname,
                cs_equation_type_t     eqtype,
                int                    dim,
                cs_param_bc_type_t     default_bc)
{
  /* Sanity checks */
  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No variable name associated to an equation structure.\n"
                " Check your initialization."));
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No equation name associated to an equation structure.\n"
                " Check your initialization."));
  if (cs_equation_by_name(eqname) != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop adding a new equation.\n"
                " Equation name %s is already defined."), eqname);

  cs_equation_t  *eq = NULL;

  BFT_MALLOC(eq, 1, cs_equation_t);

  int  eq_id = _n_equations;
  _n_equations++;
  BFT_REALLOC(_equations, _n_equations, cs_equation_t *);
  _equations[eq_id] = eq;

  switch (eqtype) {

  case CS_EQUATION_TYPE_USER:
    _n_user_equations++;
    break;

  case CS_EQUATION_TYPE_PREDEFINED:
  case CS_EQUATION_TYPE_GROUNDWATER:
    _n_predef_equations++;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This type of equation is not handled.\n"
              " Stop adding a new equation.");
    break;

  }

  eq->id = eq_id;

  /* Store eqname */
  int  len = strlen(eqname)+1;
  BFT_MALLOC(eq->name, len, char);
  strncpy(eq->name, eqname, len);

  /* Store varname */
  len = strlen(varname)+1;
  BFT_MALLOC(eq->varname, len, char);
  strncpy(eq->varname, varname, len);

  eq->param = cs_equation_create_param(eqtype, dim, default_bc);

  eq->field_id = -1;    // field is created in a second step
  eq->do_build = true;  // Force the construction of the algebraic system

  /* Set timer statistic structure to a default value */
  eq->main_ts_id = eq->solve_ts_id = -1;

  /* Algebraic system: allocated later */
  eq->matrix = NULL;
  eq->rhs = NULL;
  eq->rset = NULL;
  eq->n_sles_gather_elts = eq->n_sles_scatter_elts = 0;

  /* Builder structure for this equation */
  eq->builder = NULL;
  eq->scheme_context = NULL;

  /* Pointers of function */
  eq->init_context = NULL;
  eq->free_context = NULL;
  eq->initialize_system = NULL;
  eq->build_system = NULL;
  eq->update_field = NULL;
  eq->compute_source = NULL;
  eq->compute_flux_across_plane = NULL;
  eq->compute_cellwise_diff_flux = NULL;
  eq->postprocess = NULL;
  eq->get_extra_values = NULL;

  return  eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new user equation structure and set a first set of parameters
 *
 * \param[in] eqname        name of the equation
 * \param[in] varname       name of the variable associated to this equation
 * \param[in] dim           dimension of the unknow attached to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add_user(const char            *eqname,
                     const char            *varname,
                     int                    dim,
                     cs_param_bc_type_t     default_bc)
{
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " Empty equation name.");
  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0, " Empty variable name.");

  if ((default_bc != CS_PARAM_BC_HMG_DIRICHLET) &&
      (default_bc != CS_PARAM_BC_HMG_NEUMANN))
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of boundary condition by default.\n"
                " Valid choices are CS_PARAM_BC_HMG_DIRICHLET or"
                " CS_PARAM_BC_HMG_NEUMANN"));

  /* Add a new user equation */
  cs_equation_t  *eq =
    cs_equation_add(eqname,                // equation name
                    varname,               // variable name
                    CS_EQUATION_TYPE_USER, // type of equation
                    dim,                   // dimension of the variable
                    default_bc);           // default BC

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_destroy_all(void)
{
  if (_n_equations == 0)
    return;

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    eq->param = cs_equation_free_param(eq->param);

    /* Sanity check */
    assert(eq->matrix == NULL && eq->rhs == NULL);
    /* Since eq->rset is only shared, no free is done at this stage */

    /* Free the associated builder structure */
    cs_equation_free_builder(&(eq->builder));
    eq->scheme_context = eq->free_context(eq->scheme_context);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

    BFT_FREE(eq->name);
    BFT_FREE(eq->varname);
    BFT_FREE(eq);

  } // Loop on equations

  BFT_FREE(_equations);

  _n_equations = 0;
  _n_user_equations = 0;
  _n_predef_equations = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a synthesis of the monitoring information in the performance
 *         file
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_monitoring(void)
{
  cs_log_printf(CS_LOG_PERFORMANCE,
                "%-36s %9s %9s %9s %9s %9s %9s\n",
                " ", "SysBuild", "Diffusion", "Advection", "Reaction",
                "Source", "Extra");

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];

    /* Display high-level timer counter related to the current equation
       before deleting the structure */
    cs_equation_write_monitoring(eq->name, eq->builder);

  } // Loop on equations
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSettings for equations\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_cdo_equations          %d\n",
                _n_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_predefined_equations   %d\n",
                _n_predef_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_user_equations         %d\n",
                _n_user_equations);

  for (int  eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
    cs_log_printf(CS_LOG_SETUP,
                  "\tSummary of settings for %s eq. (variable %s)\n",
                  eq->name, eq->varname);
    cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

    cs_equation_summary_param(eq->name, eq->param);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create timer statistics structures to enable a "home-made" profiling
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_timer_stats(cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  cs_equation_param_t  *eqp = eq->param;

  /* Set timer statistics */
  if (eqp->verbosity > 0) {

    eq->main_ts_id = cs_timer_stats_create(NULL, // new root
                                           eq->name,
                                           eq->name);

    cs_timer_stats_start(eq->main_ts_id);

    if (eqp->verbosity > 1) {

      char *label = NULL;

      int  len = strlen("_solve") + strlen(eq->name) + 1;
      BFT_MALLOC(label, len, char);
      sprintf(label, "%s_solve", eq->name);
      eq->solve_ts_id = cs_timer_stats_create(eq->name, label, label);

      BFT_FREE(label);

    } // verbosity > 1

  } // verbosity > 0

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 * \param[in]  do_profiling   true or false
 *
 * \return true if all equations are steady-state otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_finalize_setup(const cs_cdo_connect_t   *connect,
                           bool                      do_profiling)
{
  if (_n_equations == 0)
    return true;

  bool  all_are_steady = true;

  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    cs_equation_param_t  *eqp = eq->param;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    if (eqp->flag & CS_EQUATION_UNSTEADY)
      all_are_steady = false;

    if (do_profiling)
      cs_equation_set_timer_stats(eq);

    /* Set function pointers */
    switch(eqp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      if (eqp->dim == 1) {

        eq->init_context = cs_cdovb_scaleq_init_context;
        eq->free_context = cs_cdovb_scaleq_free_context;
        eq->initialize_system = cs_cdovb_scaleq_initialize_system;
        eq->build_system = cs_cdovb_scaleq_build_system;
        eq->prepare_solving = _prepare_vb_solving;
        eq->update_field = cs_cdovb_scaleq_update_field;
        eq->compute_source = cs_cdovb_scaleq_compute_source;
        eq->compute_flux_across_plane =
          cs_cdovb_scaleq_compute_flux_across_plane;
        eq->compute_cellwise_diff_flux = cs_cdovb_scaleq_cellwise_diff_flux;
        eq->postprocess = cs_cdovb_scaleq_extra_op;
        eq->get_extra_values = NULL;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_vertices;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Only the scalar-valued case is handled for CDO"
                  " vertex-based schemes.\n", __func__);

      break;

    case CS_SPACE_SCHEME_CDOVCB:
      if (eqp->dim == 1) {

        eq->init_context = cs_cdovcb_scaleq_init_context;
        eq->free_context = cs_cdovcb_scaleq_free_context;
        eq->initialize_system = cs_cdovcb_scaleq_initialize_system;
        eq->build_system = cs_cdovcb_scaleq_build_system;
        eq->prepare_solving = _prepare_vb_solving;
        eq->update_field = cs_cdovcb_scaleq_update_field;
        eq->compute_source = cs_cdovcb_scaleq_compute_source;
        eq->compute_flux_across_plane =
          cs_cdovcb_scaleq_compute_flux_across_plane;
        eq->compute_cellwise_diff_flux = cs_cdovcb_scaleq_cellwise_diff_flux;
        eq->postprocess = cs_cdovcb_scaleq_extra_op;
        eq->get_extra_values = cs_cdovcb_scaleq_get_cell_values;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_vertices;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Only the scalar-valued case is handled for CDO"
                  " vertex+cell-based schemes.\n", __func__);

      break;

    case CS_SPACE_SCHEME_CDOFB:
      if (eqp->dim == 1) {

        eq->init_context = cs_cdofb_scaleq_init_context;
        eq->free_context = cs_cdofb_scaleq_free_context;
        eq->initialize_system = cs_cdofb_scaleq_initialize_system;
        eq->build_system = cs_cdofb_scaleq_build_system;
        eq->prepare_solving = _prepare_fb_solving;
        eq->update_field = cs_cdofb_scaleq_update_field;
        eq->compute_source = cs_cdofb_scaleq_compute_source;
        eq->compute_flux_across_plane = NULL;
        eq->compute_cellwise_diff_flux = NULL;
        eq->postprocess = cs_cdofb_scaleq_extra_op;
        eq->get_extra_values = cs_cdofb_scaleq_get_face_values;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_faces[0];

      }
      else if (eqp->dim == 3) {

        eq->init_context = cs_cdofb_vecteq_init_context;
        eq->free_context = cs_cdofb_vecteq_free_context;
        eq->initialize_system = cs_cdofb_vecteq_initialize_system;
        eq->build_system = cs_cdofb_vecteq_build_system;
        eq->prepare_solving = _prepare_fb_solving;
        eq->update_field = cs_cdofb_vecteq_update_field;
        eq->compute_source = cs_cdofb_vecteq_compute_source;
        eq->compute_flux_across_plane = NULL;
        eq->compute_cellwise_diff_flux = NULL;
        eq->postprocess = cs_cdofb_vecteq_extra_op;
        eq->get_extra_values = cs_cdofb_vecteq_get_face_values;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = eqp->dim * connect->n_faces[0];
        eq->n_sles_scatter_elts = eqp->dim * connect->n_faces[0];

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Only the scalar-valued and vector-valued cases are"
                  "  handled for CDO face-based schemes.\n", __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P0:
      if (eqp->dim == 1) {

        eq->init_context = cs_hho_scaleq_init_context;
        eq->free_context = cs_hho_scaleq_free_context;
        eq->initialize_system = cs_hho_scaleq_initialize_system;
        eq->build_system = cs_hho_scaleq_build_system;
        eq->prepare_solving = _prepare_fb_solving;
        eq->update_field = cs_hho_scaleq_update_field;
        eq->compute_source = cs_hho_scaleq_compute_source;
        eq->compute_flux_across_plane = NULL;
        eq->compute_cellwise_diff_flux = NULL;
        eq->postprocess = cs_hho_scaleq_extra_op;
        eq->get_extra_values = cs_hho_scaleq_get_face_values;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_faces[0];

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Only the scalar-valued case is handled for CDO"
                  " HHO schemes.\n", __func__);
      break;

    case CS_SPACE_SCHEME_HHO_P1:
      if (eqp->dim == 1) {

        eq->init_context = cs_hho_scaleq_init_context;
        eq->free_context = cs_hho_scaleq_free_context;
        eq->initialize_system = cs_hho_scaleq_initialize_system;
        eq->build_system = cs_hho_scaleq_build_system;
        eq->prepare_solving = _prepare_fb_solving;
        eq->update_field = cs_hho_scaleq_update_field;
        eq->compute_source = cs_hho_scaleq_compute_source;
        eq->compute_flux_across_plane = NULL;
        eq->compute_cellwise_diff_flux = NULL;
        eq->postprocess = cs_hho_scaleq_extra_op;
        eq->get_extra_values = cs_hho_scaleq_get_face_values;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP1];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = CS_N_FACE_DOFS_1ST * connect->n_faces[0];
        eq->n_sles_scatter_elts = CS_N_FACE_DOFS_1ST * connect->n_faces[0];

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Only the scalar-valued case is handled for CDO"
                  " HHO schemes.\n", __func__);
      break;

    case CS_SPACE_SCHEME_HHO_P2:
      if (eqp->dim == 1) {

        eq->init_context = cs_hho_scaleq_init_context;
        eq->free_context = cs_hho_scaleq_free_context;
        eq->initialize_system = cs_hho_scaleq_initialize_system;
        eq->build_system = cs_hho_scaleq_build_system;
        eq->prepare_solving = _prepare_fb_solving;
        eq->update_field = cs_hho_scaleq_update_field;
        eq->compute_source = cs_hho_scaleq_compute_source;
        eq->compute_flux_across_plane = NULL;
        eq->compute_cellwise_diff_flux = NULL;
        eq->postprocess = cs_hho_scaleq_extra_op;
        eq->get_extra_values = cs_hho_scaleq_get_face_values;

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP2];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = CS_N_FACE_DOFS_2ND * connect->n_faces[0];
        eq->n_sles_scatter_elts = CS_N_FACE_DOFS_2ND * connect->n_faces[0];

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Only the scalar-valued case is handled for CDO"
                  " HHO schemes.\n", __func__);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid scheme for the space discretization.\n"
                  " Please check your settings."));
      break;
    }

    if (cs_glob_n_ranks > 1)
      eq->n_sles_gather_elts = eq->rset->n_elts[0];

    /* Initialize cs_sles_t structure */
    cs_equation_param_set_sles(eq->name, eqp, eq->field_id);

    /* Flag this equation such that parametrization is not modifiable anymore */
    eqp->flag |= CS_EQUATION_LOCKED;

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

  return all_are_steady;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_create_fields(void)
{
  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    int  location_id = -1; // initialize values to avoid a warning
    int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

    cs_equation_t  *eq = _equations[eq_id];

    /* Sanity check */
    assert(eq != NULL);

    const cs_equation_param_t  *eqp = eq->param;

    _Bool has_previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
    if (!has_previous)
      field_mask |= CS_FIELD_STEADY;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Associate a predefined mesh_location_id to this field */
    switch (eqp->space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
    case CS_SPACE_SCHEME_CDOVCB:
      location_id = cs_mesh_location_get_id_by_name("vertices");
      break;
    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
    case CS_SPACE_SCHEME_HHO_P1:
    case CS_SPACE_SCHEME_HHO_P2:
      location_id = cs_mesh_location_get_id_by_name("cells");
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Space scheme for eq. %s is incompatible with a field.\n"
                  " Stop adding a cs_field_t structure.\n"), eq->name);
      break;
    }

    if (location_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid mesh location id (= -1) for the current field\n"));

    cs_field_t  *fld = cs_field_find_or_create(eq->varname,
                                               field_mask,
                                               location_id,
                                               eqp->dim,
                                               has_previous);

    /* Set default value for default keys */
    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);

    /* Store the related field id */
    eq->field_id = cs_field_id_by_name(eq->varname);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the builder of the algebraic system.
 *         Set the initialize condition to all variable fields associated to
 *         each cs_equation_t structure.
 *         Compute the initial source term.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_initialize(const cs_mesh_t             *mesh,
                       const cs_cdo_connect_t      *connect,
                       const cs_cdo_quantities_t   *quant,
                       const cs_time_step_t        *ts)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t *eq = _equations[i];
    assert(eq != NULL); // Sanity check

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    const cs_equation_param_t  *eqp = eq->param;

    /* Allocate and initialize a system builder */
    eq->builder = cs_equation_init_builder(eqp, mesh);
    eq->scheme_context = eq->init_context(eqp, eq->builder);

    // By default, 0 is set as initial condition
    if (eqp->n_ic_defs > 0 && ts->nt_cur < 1)
      _initialize_field_from_ic(eq);

    if (eqp->flag & CS_EQUATION_UNSTEADY)
      /* Compute the (initial) source term */
      eq->compute_source(eqp, eq->builder, eq->scheme_context);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one has to build the linear system
 *
 * \param[in]  eq        pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_build(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return false;
  else
    return eq->do_build;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system for this equation
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in]       time_step   pointer to a time step structure
 * \param[in]       dt_cur      value of the current time step
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *mesh,
                         const cs_time_step_t       *time_step,
                         double                      dt_cur,
                         cs_equation_t              *eq)
{
  assert(eq != NULL);
  CS_UNUSED(time_step);

  const cs_field_t  *fld = cs_field_by_id(eq->field_id);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Sanity checks */
  assert(eq->matrix == NULL && eq->rhs == NULL);

  /* Initialize the algebraic system to build */
  eq->initialize_system(eq->param,
                        eq->builder,
                        eq->scheme_context,
                        &(eq->matrix),
                        &(eq->rhs));

  /* Build the algebraic system to solve */
  eq->build_system(mesh, fld->val, dt_cur,
                   eq->param,
                   eq->builder,
                   eq->scheme_context,
                   eq->rhs,
                   eq->matrix);

  eq->do_build = false;

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(cs_equation_t   *eq)
{
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);
  cs_field_t  *fld = cs_field_by_id(eq->field_id);
  cs_real_t  *x = NULL, *b = NULL;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);
  if (eq->solve_ts_id > -1)
    cs_timer_stats_start(eq->solve_ts_id);

  const cs_equation_param_t  *eqp = eq->param;
  const double  r_norm = 1.0; // No renormalization by default (TODO)
  const cs_param_itsol_t  itsol_info = eqp->itsol_info;

  /* Sanity checks (up to now, only scalar field are handled) */
  assert(eq->n_sles_gather_elts <= eq->n_sles_scatter_elts);
  assert(eq->n_sles_gather_elts == cs_matrix_get_n_rows(eq->matrix));

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                " n_sles_gather_elts:  %d\n"
                " n_sles_scatter_elts: %d\n"
                " n_matrix_rows:       %d\n"
                " n_matrix_columns:    %d\n",
                eq->n_sles_gather_elts, eq->n_sles_scatter_elts,
                cs_matrix_get_n_rows(eq->matrix),
                cs_matrix_get_n_columns(eq->matrix));
#endif

  /* Handle parallelism */
  eq->prepare_solving(eq, &x, &b);

  cs_sles_convergence_state_t code = cs_sles_solve(sles,
                                                   eq->matrix,
                                                   CS_HALO_ROTATION_IGNORE,
                                                   itsol_info.eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   b,
                                                   x,
                                                   0,      // aux. size
                                                   NULL);  // aux. buffers

  if (eq->param->sles_verbosity > 0) {

    const cs_lnum_t  size = eq->n_sles_gather_elts;
    const cs_lnum_t  *row_index, *col_id;
    const cs_real_t  *d_val, *x_val;

    cs_matrix_get_msr_arrays(eq->matrix, &row_index, &col_id, &d_val, &x_val);

    cs_gnum_t  nnz = row_index[size];
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&nnz, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  "  <%s/sles_cvg> code  %d n_iters  %d residual  % -8.4e"
                  " nnz %lu\n",
                  eq->name, code, n_iters, residual, nnz);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_DBG > 1
    cs_dbg_darray_to_listing("EQ %s: AFTER.SOLVE >> X", eq->name, size, x, 9);
    cs_dbg_darray_to_listing("EQ %s: SOLVE >> RHS", eq->name, size, b, 9);
#if CS_EQUATION_DBG > 2
    cs_dbg_iarray_to_listing("EQ %s: ROW_INDEX",
                             eq->name, size + 1, row_index, 9);
    cs_dbg_iarray_to_listing("EQ %s: COLUMN_ID", eq->name, nnz, col_id, 9);
    cs_dbg_darray_to_listing("EQ %s: D_VAL", eq->name, size, d_val, 9);
    cs_dbg_darray_to_listing("EQ %s: X_VAL", eq->name, nnz, x_val, 9);
#endif
#endif
  }

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_range_set_scatter(eq->rset,
                         CS_REAL_TYPE, 1, // type and stride
                         x,
                         x);

    cs_range_set_scatter(eq->rset,
                         CS_REAL_TYPE, 1, // type and stride
                         b,
                         eq->rhs);

  }

  if (eq->solve_ts_id > -1)
    cs_timer_stats_stop(eq->solve_ts_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_DBG > 1
  cs_dbg_array_fprintf(NULL, "sol.log", 1e-16, eq->n_sles_gather_elts, x, 6);
  cs_dbg_array_fprintf(NULL, "rhs.log", 1e-16, eq->n_sles_gather_elts, b, 6);
#endif

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Define the new field value for the current time */
  eq->update_field(x, eq->rhs, eq->param,
                   eq->builder, eq->scheme_context, fld->val);

  if (eq->param->flag & CS_EQUATION_UNSTEADY)
    eq->do_build = true; /* Improvement: exhibit cases where a new build
                            is not needed */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  /* Free memory */
  BFT_FREE(x);
  if (b != eq->rhs)
    BFT_FREE(b);
  BFT_FREE(eq->rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&(eq->matrix));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each face of the mesh for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the face values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  if (eq->param->space_scheme == CS_SPACE_SCHEME_CDOFB)
    return eq->get_extra_values(eq->scheme_context);

  else if (eq->param->space_scheme == CS_SPACE_SCHEME_HHO_P0 ||
           eq->param->space_scheme == CS_SPACE_SCHEME_HHO_P1 ||
           eq->param->space_scheme == CS_SPACE_SCHEME_HHO_P2) {

    if (eq->param->dim == 1)    /* scalar-valued unknown */
      return cs_hho_scaleq_get_face_values(eq->scheme_context);
    else
      return NULL;               /* TODO */

  }
  else {

    if (eq->get_extra_values == NULL) {
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: No function defined for this operation in eq. %s"),
                __func__, eq->name);
      return NULL;                /* TODO */

    }

  }

  return NULL; /* Avoid a warning */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each cell centers for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the cell values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  switch (eq->param->space_scheme) {
  case CS_SPACE_SCHEME_CDOFB:
    {
      cs_field_t  *fld = cs_field_by_id(eq->field_id);

      return fld->val;
    }
    break;

  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    {
      if (eq->param->dim == 1)    /* scalar-valued unknown */
        return cs_hho_scaleq_get_cell_values(eq->scheme_context);
      else
        return NULL;               /* TODO */
    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    return eq->get_extra_values(eq->scheme_context);
    break;

  default:
    {
      if (eq->get_extra_values == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: No function defined for this operation in eq. %s"),
                  __func__, eq->name);
      return NULL; // Avoid aa warning
    }
    break;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux accross a plane defined
 *         by a mesh location structure attached to the name ml_name.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      ml_name     name of the related mesh location
 * \param[in]      direction   vector indicating in which direction flux is > 0
 * \param[in, out] diff_flux   value of the diffusive part of the flux
 * \param[in, out] conv_flux   value of the convective part of the flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_flux_across_plane(const cs_equation_t   *eq,
                                      const char            *ml_name,
                                      const cs_real_3_t      direction,
                                      cs_real_t             *diff_flux,
                                      cs_real_t             *conv_flux)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->compute_flux_across_plane == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _(" Computation of the diffusive and convective flux across\n"
                " a plane is not available for equation %s\n"), eq->name);
    return; // Avoid a warning
  }

  /* Get the mesh location id from its name */
  const int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Perform the computation */
  eq->compute_flux_across_plane(direction,
                                fld->val,
                                ml_id,
                                eq->param,
                                eq->builder,
                                eq->scheme_context,
                                diff_flux, conv_flux);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all cell faces.
 *         Primal or dual faces are considered according to the space scheme.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      location    indicate where the flux has to be computed
 * \param[in, out] diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diff_flux_cellwise(const cs_equation_t   *eq,
                                       cs_flag_t              location,
                                       cs_real_t             *diff_flux)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->compute_cellwise_diff_flux == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _(" Cellwise computation of the diffusive flux is not\n"
                " available for equation %s\n"), eq->name);
    return; // Avoid a warning
  }

  if (eq->builder == NULL)
    return;

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Perform the computation */
  eq->compute_cellwise_diff_flux(fld->val,
                                 eq->param,
                                 eq->builder,
                                 eq->scheme_context,
                                 location,
                                 diff_flux);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the discrete gradient at vertices
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in, out] v_gradient  gradient at vertices
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_vtx_field_gradient(const cs_equation_t   *eq,
                                       cs_real_t             *v_gradient)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);
  assert(v_gradient != NULL);

  const cs_equation_param_t  *eqp = eq->param;

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVCB:
    cs_cdovcb_scaleq_vtx_gradient(fld->val, eq->builder, eq->scheme_context,
                                  v_gradient);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of scheme for compting the gradient at vertices");
    break;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to all equations
 *
 * \param[in]  ts      pointer to a cs_time_step_t struct.
 * \param[in]  dt      value of the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_post_all(const cs_time_step_t    *ts,
                           double                   dt)
{
  if (_n_equations < 1)
    return;

  CS_UNUSED(dt);

  int  len;

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    char *postlabel = NULL;

    const cs_field_t  *field = cs_field_by_id(eq->field_id);
    const cs_equation_param_t  *eqp = eq->param;

    /* Cases where a post-processing is not required */
    if (eqp->process_flag == 0)
      continue;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Post-processing of a common adimensionnal quantities: the Peclet
       number */
    if (eqp->process_flag & CS_EQUATION_POST_PECLET) {

      len = strlen(eq->name) + 7 + 1;
      BFT_MALLOC(postlabel, len, char);
      sprintf(postlabel, "%s.Peclet", eq->name);

      /* Compute the Peclet number in each cell */
      double  *peclet = cs_equation_get_tmpbuf();
      cs_advection_get_peclet(eqp->adv_field, eqp->diffusion_property, peclet);

      /* Post-process */
      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        postlabel,
                        1,
                        true,           // interlace
                        true,           // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        peclet,         // values on cells
                        NULL,           // values at internal faces
                        NULL,           // values at border faces
                        ts);            // time step management struct.

      BFT_FREE(postlabel);

    } // Peclet number

    /* Perform post-processing specific to a numerical scheme */
    eq->postprocess(eq->name, field,
                    eq->param,
                    eq->builder,
                    eq->scheme_context);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
