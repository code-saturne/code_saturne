/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_cdovb_vecteq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdoeb_vecteq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_field_default.h"
#include "cs_hho_scaleq.h"
#include "cs_hho_vecteq.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_range_set.h"
#include "cs_sles.h"
#include "cs_timer_stats.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

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
  N_(" %s: Stop setting an empty cs_equation_t structure.\n"
     " Please check your settings.\n");

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values for the variable related to an equation
 *
 * \param[in]       eq       pointer to a cs_equation_t structure
 * \param[in]       ts       pointer to cs_time_step_t structure
 * \param[in]       tag      tag to add to the equation name to build the label
 * \param[in, out]  label    label for the postprocessing
 * \param[in]       values   pointer to the array of values to post
 */
/*----------------------------------------------------------------------------*/

static inline void
_post_balance_at_vertices(const cs_equation_t   *eq,
                          const cs_time_step_t  *ts,
                          const char            *tag,
                          char                  *label,
                          const cs_real_t       *values)

{
  sprintf(label, "%s.Balance.%s", eq->param->name, tag);

  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           eq->param->dim,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           values,
                           ts);
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

  const cs_real_t  *f_values = eq->get_face_values(eq->scheme_context);
  const int  stride = 1;  /* Since the global numbering is adapted in each
                             case (scalar-, vector-valued equations) */

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
                        CS_REAL_TYPE,  /* type */
                        stride,        /* stride */
                        f_values,      /* in: size = n_sles_scatter_elts */
                        x);            /* out: size = n_sles_gather_elts */

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
                         eq->n_sles_scatter_elts, stride, false, CS_REAL_TYPE,
                         b);

    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE,  /* type */
                        stride,        /* stride */
                        b,             /* in: size = n_sles_scatter_elts */
                        b);            /* out: size = n_sles_gather_elts */

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the pointers of function for the given equation.
 *         Case of scalar-valued HHO schemes
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_scal_hho_function_pointers(cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  eq->init_context = cs_hho_scaleq_init_context;
  eq->free_context = cs_hho_scaleq_free_context;

  /* New functions */
  eq->init_field_values = cs_hho_scaleq_init_values;
  eq->solve = NULL;
  eq->solve_steady_state = NULL;

  eq->postprocess = cs_hho_scaleq_extra_op;
  eq->read_restart = cs_hho_scaleq_read_restart;
  eq->write_restart = cs_hho_scaleq_write_restart;

  /* Function pointers to retrieve values at mesh locations */
  eq->get_vertex_values = NULL;
  eq->get_cell_values = cs_hho_scaleq_get_cell_values;
  eq->get_face_values = cs_hho_scaleq_get_face_values;

  /* Deprecated functions */
  eq->initialize_system = cs_hho_scaleq_initialize_system;
  eq->set_dir_bc = NULL;
  eq->build_system = cs_hho_scaleq_build_system;
  eq->prepare_solving = _prepare_fb_solving;
  eq->update_field = cs_hho_scaleq_update_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the pointers of function for the given equation.
 *         Case of vector-valued HHO schemes
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_vect_hho_function_pointers(cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  eq->init_context = cs_hho_vecteq_init_context;
  eq->free_context = cs_hho_vecteq_free_context;

  /* New functions */
  eq->init_field_values = cs_hho_vecteq_init_values;
  eq->solve = NULL;
  eq->solve_steady_state = NULL;

  eq->postprocess = cs_hho_vecteq_extra_op;
  eq->read_restart = cs_hho_vecteq_read_restart;
  eq->write_restart = cs_hho_vecteq_write_restart;

  /* Function pointers to retrieve values at mesh locations */
  eq->get_vertex_values = NULL;
  eq->get_cell_values = cs_hho_vecteq_get_cell_values;
  eq->get_face_values = cs_hho_vecteq_get_face_values;

  /* Deprecated functions */
  eq->initialize_system = cs_hho_vecteq_initialize_system;
  eq->build_system = cs_hho_vecteq_build_system;
  eq->prepare_solving = _prepare_fb_solving;
  eq->update_field = cs_hho_vecteq_update_field;
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
    assert(_eq != NULL);
    cs_equation_param_t  *eqp = _eq->param;
    assert(eqp != NULL);
    if (strlen(eqp->name) == len_in)
      if (strcmp(eqname, eqp->name) == 0)
        return _eq;

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the asociated field to a \ref cs_equation_t structure
 *         has name equal to fld_name
 *
 * \param[in]  eq          pointer to a \ref cs_equation_t structure to test
 * \param[in]  fld_name    name of the field
 *
 * \return true if the \ref cs_equation_t structure has an associated field
 *         named fld_name, otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_has_field_name(const cs_equation_t  *eq,
                           const char           *fld_name)
{
  if (eq == NULL)
    return false;

  cs_field_t  *fld = cs_field_by_id(eq->field_id);
  if (fld == NULL)
    return false;

  if (strcmp(fld->name, fld_name) == 0)
    return true;
  else
    return false;
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
 * \brief  Find the \ref cs_equation_t structure with id eq_id
 *         Return NULL if not find
 *
 * \param[in]  eq_id    id of the equation to find
 *
 * \return a pointer to a \ref cs_equation_t structure or NULL if not found
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
  cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  if (eqp == NULL)
    return NULL;
  else
    return eqp->name;
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
 * \brief  Return the id related to the variable field structure associated to
 *         the cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return an integer (-1 if the field is not defined)
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_field_id(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return -1;
  else
    return eq->field_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the field structure for the (normal) boundary flux associated
 *         to a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_field_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_equation_get_boundary_flux(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else if (eq->boundary_flux_id < 0)
    return NULL;
  else
    return cs_field_by_id(eq->boundary_flux_id);
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
 * \brief  Redefine the flag associated to an equation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 * \param[in]       flag     new flag to set
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_set_flag(cs_equation_t    *eq,
                     cs_flag_t         flag)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq, __func__);
  assert(eq->param != NULL);

  eq->param->flag = flag;;
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
 *         reaction term with id equal to reaction_id and related to this
 *         equation
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
 * \brief  Return the type of numerical scheme used for the discretization in
 *         time
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_param_time_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_param_time_scheme_t
cs_equation_get_time_scheme(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_TIME_N_SCHEMES;
  else if (eq->param == NULL)
    return CS_TIME_N_SCHEMES;
  else
    return eq->param->time_scheme;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the value of the theta parameter in theta time scheme
 *         discretization
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the value of the theta coefficient. -1 if the time scheme is not
 *          a theta time scheme
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_equation_get_theta_time_val(const cs_equation_t    *eq)
{
  cs_real_t  theta = -1;

  if (eq == NULL)
    return theta;
  else if (eq->param == NULL)
    return theta;

  else {

    switch (eq->param->time_scheme) {

    case CS_TIME_SCHEME_THETA:
      theta = eq->param->theta;
      break;
    case CS_TIME_SCHEME_CRANKNICO:
      theta = 0.5;
      break;
    case CS_TIME_SCHEME_EULER_IMPLICIT:
      theta = 1;
      break;
    case CS_TIME_SCHEME_EULER_EXPLICIT:
      theta = 0.;
      break;

    default:
      break;
    }

  }

  return theta;
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
 * \brief  Return true is the given equation is steady otherwise false
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_uses_new_mechanism(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return false;
  assert(eq->param != NULL);

  if (eq->param->dim == 1) {
    if ((eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB)  ||
        (eq->param->space_scheme == CS_SPACE_SCHEME_CDOVCB) ||
        (eq->param->space_scheme == CS_SPACE_SCHEME_CDOFB))
      return true;
  }
  else if (eq->param->dim == 3) {
    if ((eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB) ||
        (eq->param->space_scheme == CS_SPACE_SCHEME_CDOFB))
      return true;
  }

  return false;
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
              _(" %s: No variable name associated to an equation structure.\n"
                " Check your initialization."), __func__);
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s No equation name associated to an equation structure.\n"
                " Check your initialization."), __func__);
  if (cs_equation_by_name(eqname) != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Stop adding a new equation.\n"
                " Equation name %s is already defined."), __func__, eqname);

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
  case CS_EQUATION_TYPE_MAXWELL:
  case CS_EQUATION_TYPE_NAVSTO:
    _n_predef_equations++;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: This type of equation is not handled.\n"
              " Stop adding a new equation.", __func__);
    break;

  }

  eq->id = eq_id;

  /* Store varname */
  int  len = strlen(varname)+1;
  BFT_MALLOC(eq->varname, len, char);
  strncpy(eq->varname, varname, len);

  eq->param = cs_equation_create_param(eqname, eqtype, dim, default_bc);

  eq->field_id = -1;           /* This field is created in a second step */
  eq->boundary_flux_id = -1;   /* Not always defined (done in a second step) */

  /* Algebraic system: allocated later */
  eq->rset = NULL;
  eq->n_sles_gather_elts = eq->n_sles_scatter_elts = 0;

  /* Builder structure for this equation */
  eq->builder = NULL;
  eq->scheme_context = NULL;

  /* Pointers of function */
  eq->init_context = NULL;
  eq->free_context = NULL;

  /* Postprocessing */
  eq->compute_balance = NULL;
  eq->postprocess = NULL;

  /* Restart */
  eq->read_restart = NULL;
  eq->write_restart = NULL;

  /* Function pointers to retrieve values at mesh locations */
  eq->get_vertex_values = NULL;
  eq->get_cell_values = NULL;
  eq->get_face_values = NULL;

  /* New functions */
  eq->init_field_values = NULL;
  eq->solve = NULL;
  eq->solve_steady_state = NULL;

  /* Deprecated functions */
  eq->initialize_system = NULL;
  eq->set_dir_bc = NULL;
  eq->build_system = NULL;
  eq->update_field = NULL;

  /* Deprecated members related to deprecated functions */
  eq->matrix = NULL;
  eq->rhs = NULL;

  /* Set timer statistic structure to a default value */
  eq->main_ts_id = cs_timer_stats_create(NULL, /* new root */
                                         eqname,
                                         eqname);

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
    bft_error(__FILE__, __LINE__, 0, " %s: Empty equation name.", __func__);
  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Empty variable name.", __func__);

  if ((default_bc != CS_PARAM_BC_HMG_DIRICHLET) &&
      (default_bc != CS_PARAM_BC_HMG_NEUMANN))
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid type of boundary condition by default.\n"
                " Valid choices are CS_PARAM_BC_HMG_DIRICHLET or"
                " CS_PARAM_BC_HMG_NEUMANN"), __func__);

  /* Add a new user equation */
  cs_equation_t  *eq =
    cs_equation_add(eqname,                /* equation name */
                    varname,               /* variable name */
                    CS_EQUATION_TYPE_USER, /* type of equation */
                    dim,                   /* dimension of the variable */
                    default_bc);           /* default BC */

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

    BFT_FREE(eq->varname);
    BFT_FREE(eq);

  } /* Loop on equations */

  BFT_FREE(_equations);

  _n_equations = 0;
  _n_user_equations = 0;
  _n_predef_equations = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a steady-state computation is requested according to the
 *         setting
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_steady_state_solve(void)
{
  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];

    if (cs_equation_is_steady(eq))
      return true;

  } /* Loop on equations */

  return false;
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
    cs_equation_write_monitoring(eq->param->name, eq->builder);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the count of equations of each macro type
 *
 * \param[out]  n_equations          total number of equations
 * \param[out]  n_predef_equations   number of predefined equations
 * \param[out]  n_user_equations     number of user equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_get_count(int      *n_equations,
                      int      *n_predef_equations,
                      int      *n_user_equations)
{
  *n_equations = _n_equations;
  *n_predef_equations = _n_predef_equations;
  *n_user_equations = _n_user_equations;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP, "\nSettings for equations\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", h1_sep);

  for (int  eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    assert(eq->param != NULL);

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    cs_log_printf(CS_LOG_SETUP,
                  "\nSummary of settings for %s eq. (variable %s)\n",
                  eq->param->name, eq->varname);
    cs_log_printf(CS_LOG_SETUP, "%s", h2_sep);

    cs_equation_summary_param(eq->param);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname for the default settigns
 *
 * \param[in] key      key related to the member of eq to set
 * \param[in] keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_default_param(cs_equation_key_t      key,
                              const char            *keyval)
{
  if (_n_equations == 0)
    return;

  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    if (eq == NULL)
      continue;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    cs_equation_set_param(eq->param, key, keyval);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the linear algebra requirements
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_sles(void)
{
  if (_n_equations == 0)
    return;

  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    cs_equation_param_t  *eqp = eq->param;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Initialize cs_sles_t structure */
    if (eqp->type != CS_EQUATION_TYPE_NAVSTO)
      cs_equation_param_set_sles(eqp);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared structures among the activated class of discretization
 *         schemes
 *
 * \param[in]  connect          pointer to a cs_cdo_connect_t structure
 * \param[in]  quant            pointer to additional mesh quantities struct.
 * \param[in]  time_step        pointer to a time step structure
 * \param[in]  eb_scheme_flag   metadata for Eb schemes
 * \param[in]  fb_scheme_flag   metadata for Fb schemes
 * \param[in]  vb_scheme_flag   metadata for Vb schemes
 * \param[in]  vcb_scheme_flag  metadata for V+C schemes
 * \param[in]  hho_scheme_flag  metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_shared_structures(const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  const cs_time_step_t        *time_step,
                                  cs_flag_t                    eb_scheme_flag,
                                  cs_flag_t                    fb_scheme_flag,
                                  cs_flag_t                    vb_scheme_flag,
                                  cs_flag_t                    vcb_scheme_flag,
                                  cs_flag_t                    hho_scheme_flag)
{
  if (vb_scheme_flag > 0 || vcb_scheme_flag > 0) {

    if (vb_scheme_flag  & CS_FLAG_SCHEME_SCALAR ||
        vcb_scheme_flag & CS_FLAG_SCHEME_SCALAR) {

      cs_matrix_structure_t  *ms
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_VTX_SCAL);

      if (vb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
        cs_cdovb_scaleq_init_common(quant, connect, time_step, ms);

      if (vcb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
        cs_cdovcb_scaleq_init_common(quant, connect, time_step, ms);

    } /* scalar-valued DoFs */

    if (vb_scheme_flag  & CS_FLAG_SCHEME_VECTOR ||
        vcb_scheme_flag & CS_FLAG_SCHEME_VECTOR) {

      cs_matrix_structure_t  *ms
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_VTX_VECT);

      if (vb_scheme_flag & CS_FLAG_SCHEME_VECTOR)
        cs_cdovb_vecteq_init_common(quant, connect, time_step, ms);

      /* TODO: VCb schemes DoFs */

    } /* vector-valued DoFs */

  } /* Vertex-based class of discretization schemes */

  if (eb_scheme_flag > 0) {

    if (eb_scheme_flag  & CS_FLAG_SCHEME_SCALAR) {

      cs_matrix_structure_t  *ms
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_EDGE_SCAL);

      cs_cdoeb_vecteq_init_common(quant, connect, time_step, ms);

    } /* scalar-valued DoFs */

  } /* Edge-based class of discretization schemes */

  if (fb_scheme_flag > 0 || hho_scheme_flag > 0) {

    if (cs_flag_test(fb_scheme_flag,
                     CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR)) {

      cs_matrix_structure_t  *ms
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_SP0);

      cs_cdofb_scaleq_init_common(quant, connect, time_step, ms);

    } /* Scalar-valued CDO-Fb DoFs */

    if (cs_flag_test(fb_scheme_flag,
                     CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR)) {

      cs_matrix_structure_t  *ms
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_SP1);

      cs_cdofb_vecteq_init_common(quant, connect, time_step, ms);

    } /* Vector-valued CDO-Fb DoFs */

    if (hho_scheme_flag & CS_FLAG_SCHEME_SCALAR) {

      cs_matrix_structure_t  *ms0
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_SP0);
      cs_matrix_structure_t  *ms1
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_SP1);
      cs_matrix_structure_t  *ms2
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_SP2);

      cs_hho_scaleq_init_common(hho_scheme_flag,
                                quant, connect, time_step,
                                ms0, ms1, ms2);

    } /* Scalar-valued HHO schemes DoFs */

    if (hho_scheme_flag & CS_FLAG_SCHEME_VECTOR) {

      cs_matrix_structure_t  *ms0
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_VHP0);
      cs_matrix_structure_t  *ms1
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_VHP1);
      cs_matrix_structure_t  *ms2
        = cs_equation_get_matrix_structure(CS_CDO_CONNECT_FACE_VHP2);

      cs_hho_vecteq_init_common(hho_scheme_flag,
                                quant, connect, time_step,
                                ms0, ms1, ms2);

    } /* Vector-valued HHO schemes DoFs */

  } /* Face-based class of discretization schemes */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Release shared structures among the activated class of discretization
 *         schemes
 *
 * \param[in]  vb_scheme_flag   metadata for Vb schemes
 * \param[in]  vcb_scheme_flag  metadata for V+C schemes
 * \param[in]  fb_scheme_flag   metadata for Fb schemes
 * \param[in]  hho_scheme_flag  metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_unset_shared_structures(cs_flag_t    vb_scheme_flag,
                                    cs_flag_t    vcb_scheme_flag,
                                    cs_flag_t    fb_scheme_flag,
                                    cs_flag_t    hho_scheme_flag)
{
  /* Free common structures specific to a numerical scheme */
  if (vb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdovb_scaleq_finalize_common();

  if (vb_scheme_flag & CS_FLAG_SCHEME_VECTOR)
    cs_cdovb_vecteq_finalize_common();

  if (vcb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdovcb_scaleq_finalize_common();

  if (fb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdofb_scaleq_finalize_common();

  if (fb_scheme_flag & CS_FLAG_SCHEME_VECTOR)
    cs_cdofb_vecteq_finalize_common();

  if (hho_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_hho_scaleq_finalize_common();

  if (hho_scheme_flag & CS_FLAG_SCHEME_VECTOR)
    cs_hho_vecteq_finalize_common();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a \ref cs_range_set_t structures for synchronization when
 *         computing in parallel mode.
 *
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_range_set(const cs_cdo_connect_t   *connect)

{
  if (_n_equations == 0)
    return;

  const char  s_err_msg[] =
    "%s: Only the scalar-valued case is handled for this scheme.\n";
  const char  sv_err_msg[] =
    "%s: Only the scalar-valued and vector-valued case are handled"
    "for this scheme.\n";

  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];

  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    cs_equation_param_t  *eqp = eq->param;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Set function pointers */
    switch(eqp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      if (eqp->dim == 1) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_vertices;

      }
      else if (eqp->dim == 3) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_VTX_VECT];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eqp->dim * connect->n_vertices;
        eq->n_sles_scatter_elts = eqp->dim * connect->n_vertices;

      }
      else
        bft_error(__FILE__, __LINE__, 0, sv_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_CDOVCB:
      if (eqp->dim == 1) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_vertices;

      }
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_CDOFB:
      if (eqp->dim == 1) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = n_faces;

      }
      else if (eqp->dim == 3) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = eqp->dim * n_faces;
        eq->n_sles_scatter_elts = eqp->dim * n_faces;

      }
      else
        bft_error(__FILE__, __LINE__, 0, sv_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P0:
      if (eqp->dim == 1) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];

        /* Set the size of the algebraic system arising from the cellwise
           process */
        eq->n_sles_gather_elts = eq->n_sles_scatter_elts = n_faces;

      }
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P1:
      if (eqp->dim == 1) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP1];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = CS_N_FACE_DOFS_1ST * n_faces;
        eq->n_sles_scatter_elts = CS_N_FACE_DOFS_1ST * n_faces;

      }
      else if (eqp->dim == 3) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_VHP1];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = 3*CS_N_FACE_DOFS_1ST * n_faces;
        eq->n_sles_scatter_elts = 3*CS_N_FACE_DOFS_1ST * n_faces;

      }
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P2:
      if (eqp->dim == 1) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP2];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = CS_N_FACE_DOFS_2ND * n_faces;
        eq->n_sles_scatter_elts = CS_N_FACE_DOFS_2ND * n_faces;

      }
      else if (eqp->dim == 3) {

        /* Set the cs_range_set_t structure */
        eq->rset = connect->range_sets[CS_CDO_CONNECT_FACE_VHP2];

        /* Set the size of the algebraic system arising from the cellwise
           process (OK for a sequential run) */
        eq->n_sles_gather_elts = 3*CS_N_FACE_DOFS_2ND * n_faces;
        eq->n_sles_scatter_elts = 3*CS_N_FACE_DOFS_2ND * n_faces;

      }
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid scheme for the space discretization.\n"
                  " Please check your settings."));
      break;
    }

    if (cs_glob_n_ranks > 1)
      eq->n_sles_gather_elts = eq->rset->n_elts[0];

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *         After this call, parameters related to an equation are set once for
 *         all
 *
 * \return true if all equations are steady-state otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_set_functions(void)
{
  if (_n_equations == 0)
    return true;

  const char  s_err_msg[] =
    "%s: Only the scalar-valued case is handled for this scheme.\n";
  const char  sv_err_msg[] =
    "%s: Only the scalar-valued and vector-valued case are handled"
    "for this scheme.\n";

  bool  all_are_steady = true;

  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    cs_equation_param_t  *eqp = eq->param;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    if (eqp->flag & CS_EQUATION_UNSTEADY)
      all_are_steady = false;
    else
      cs_equation_set_param(eqp, CS_EQKEY_TIME_SCHEME, "steady");

    /* Apply the last modifications to the cs_equation_param_t structure */
    cs_equation_param_last_stage(eqp);

    /* Set function pointers */
    switch(eqp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      if (eqp->dim == 1) {

        eq->init_context = cs_cdovb_scaleq_init_context;
        eq->free_context = cs_cdovb_scaleq_free_context;
        eq->init_field_values = cs_cdovb_scaleq_init_values;

        /* deprecated */
        eq->set_dir_bc = NULL;
        eq->initialize_system = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */
        eq->solve_steady_state = cs_cdovb_scaleq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
          eq->solve = eq->solve_steady_state;
          break;
        case CS_TIME_SCHEME_EULER_IMPLICIT:
          eq->solve = cs_cdovb_scaleq_solve_implicit;
          break;
        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
          eq->solve = cs_cdovb_scaleq_solve_theta;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->compute_balance = cs_cdovb_scaleq_balance;
        eq->postprocess = cs_cdovb_scaleq_extra_op;

        eq->read_restart = NULL;
        eq->write_restart = NULL;

        /* Function pointers to retrieve values at mesh locations */
        eq->get_vertex_values = cs_cdovb_scaleq_get_vertex_values;
        eq->get_cell_values = cs_cdovb_scaleq_get_cell_values;
        eq->get_face_values = NULL;

        eq->get_cw_build_structures = cs_cdovb_scaleq_get;

      }
      else if (eqp->dim == 3) {

        eq->init_context = cs_cdovb_vecteq_init_context;
        eq->free_context = cs_cdovb_vecteq_free_context;
        eq->init_field_values = cs_cdovb_vecteq_init_values;

        /* Deprecated */
        eq->set_dir_bc = NULL;
        eq->initialize_system = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */
        eq->solve_steady_state = cs_cdovb_vecteq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
        eq->solve = eq->solve_steady_state;
          break;
        case CS_TIME_SCHEME_EULER_IMPLICIT:
        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdovb_vecteq_extra_op;

        eq->read_restart = NULL;
        eq->write_restart = NULL;

        /* Function pointers to retrieve values at mesh locations */
        eq->get_vertex_values = cs_cdovb_vecteq_get_vertex_values;
        eq->get_cell_values = cs_cdovb_vecteq_get_cell_values;
        eq->get_face_values = NULL;

        eq->get_cw_build_structures = cs_cdovb_vecteq_get;

      }
      else
        bft_error(__FILE__, __LINE__, 0, sv_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_CDOVCB:
      if (eqp->dim == 1) {

        eq->init_context = cs_cdovcb_scaleq_init_context;
        eq->free_context = cs_cdovcb_scaleq_free_context;
        eq->init_field_values = cs_cdovcb_scaleq_init_values;

        /* Deprecated */
        eq->set_dir_bc = NULL;
        eq->initialize_system = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */
        eq->solve_steady_state = cs_cdovcb_scaleq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
          eq->solve = eq->solve_steady_state;
          break;

        case CS_TIME_SCHEME_EULER_IMPLICIT:
          eq->solve = cs_cdovcb_scaleq_solve_implicit;
          break;

        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
          eq->solve = cs_cdovcb_scaleq_solve_theta;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdovcb_scaleq_extra_op;
        eq->read_restart = cs_cdovcb_scaleq_read_restart;
        eq->write_restart = cs_cdovcb_scaleq_write_restart;

        /* Function pointers to retrieve values at mesh locations */
        eq->get_vertex_values = cs_cdovcb_scaleq_get_vertex_values;
        eq->get_cell_values = cs_cdovcb_scaleq_get_cell_values;
        eq->get_face_values = NULL;

        eq->get_cw_build_structures = cs_cdovcb_scaleq_get;

      }
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_CDOFB:
      if (eqp->dim == 1) {

        eq->init_context = cs_cdofb_scaleq_init_context;
        eq->free_context = cs_cdofb_scaleq_free_context;
        eq->init_field_values = cs_cdofb_scaleq_init_values;

        /* Deprecated */
        eq->set_dir_bc = NULL;
        eq->initialize_system = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */
        eq->solve_steady_state = cs_cdofb_scaleq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
          eq->solve = eq->solve_steady_state;
          break;

        case CS_TIME_SCHEME_EULER_IMPLICIT:
          eq->solve = cs_cdofb_scaleq_solve_implicit;
          break;

        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
          eq->solve = cs_cdofb_scaleq_solve_theta;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->compute_balance = cs_cdofb_scaleq_balance;
        eq->postprocess = cs_cdofb_scaleq_extra_op;

        eq->read_restart = cs_cdofb_scaleq_read_restart;
        eq->write_restart = cs_cdofb_scaleq_write_restart;

        /* Function pointers to retrieve values at mesh locations */
        eq->get_vertex_values = NULL;
        eq->get_cell_values = cs_cdofb_scaleq_get_cell_values;
        eq->get_face_values = cs_cdofb_scaleq_get_face_values;

        eq->get_cw_build_structures = cs_cdofb_scaleq_get;

      }
      else if (eqp->dim == 3) {

        eq->init_context = cs_cdofb_vecteq_init_context;
        eq->free_context = cs_cdofb_vecteq_free_context;
        eq->init_field_values = cs_cdofb_vecteq_init_values;

        /* Depreacted */
        eq->initialize_system = NULL;
        eq->set_dir_bc = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */
        eq->solve_steady_state = cs_cdofb_vecteq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
          eq->solve = eq->solve_steady_state;
          break;

        case CS_TIME_SCHEME_EULER_IMPLICIT:
          eq->solve = cs_cdofb_vecteq_solve_implicit;
          break;

        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
          eq->solve = cs_cdofb_vecteq_solve_theta;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdofb_vecteq_extra_op;

        eq->read_restart = cs_cdofb_vecteq_read_restart;
        eq->write_restart = cs_cdofb_vecteq_write_restart;

        /* Function pointers to retrieve values at mesh locations */
        eq->get_vertex_values = NULL;
        eq->get_cell_values = cs_cdofb_vecteq_get_cell_values;
        eq->get_face_values = cs_cdofb_vecteq_get_face_values;

        eq->get_cw_build_structures = cs_cdofb_vecteq_get;
      }
      else
        bft_error(__FILE__, __LINE__, 0, sv_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P0:
      if (eqp->dim == 1) /* Set pointers of function */
        _set_scal_hho_function_pointers(eq);
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P1:
      if (eqp->dim == 1) /* Set pointers of function */
        _set_scal_hho_function_pointers(eq);

      else if (eqp->dim == 3) /* Set pointers of function */
        _set_vect_hho_function_pointers(eq);

      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P2:
      if (eqp->dim == 1) /* Set pointers of function */
        _set_scal_hho_function_pointers(eq);

      else if (eqp->dim == 3) /* Set pointers of function */
        _set_vect_hho_function_pointers(eq);

      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Eq. %s. Invalid scheme for the space discretization.\n"
                  " Please check your settings."), __func__, eqp->name);
      break;
    }

    /* Flag this equation such that parametrization is not modifiable anymore */
    eqp->flag |= CS_EQUATION_LOCKED;

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */

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

    int  location_id = -1; /* initialize values to avoid a warning */
    cs_equation_t  *eq = _equations[eq_id];

    /* Sanity check */
    assert(eq != NULL);

    cs_equation_param_t  *eqp = eq->param;

    /* Redondant definition to handle C/FORTRAN */
    int  has_previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? 1 : 0;
    bool  previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Associate a predefined mesh_location_id to this field */
    switch (eqp->space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
    case CS_SPACE_SCHEME_CDOVCB:
      location_id = cs_mesh_location_get_id_by_name("vertices");
      break;
    case CS_SPACE_SCHEME_CDOEB:
    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
    case CS_SPACE_SCHEME_HHO_P1:
    case CS_SPACE_SCHEME_HHO_P2:
      location_id = cs_mesh_location_get_id_by_name("cells");
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Space scheme for eq. %s is incompatible with a field.\n"
                  " Stop adding a cs_field_t structure.\n"), eqp->name);
      break;
    }

    if (location_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid mesh location id (= -1) for the current field\n"));

    /* Store the related field id */
    eq->field_id = cs_variable_cdo_field_create(eq->varname,
                                                NULL, /* label */
                                                location_id,
                                                eqp->dim,
                                                has_previous);

    /* SLES is associated to a field_id */
    eqp->sles_param.field_id = eq->field_id;

    if (eqp->process_flag & CS_EQUATION_POST_NORMAL_FLUX) {

      /* Add a field for the normal boundary flux */
      location_id = cs_mesh_location_get_id_by_name("boundary_faces");

      char  *bdy_flux_name = NULL;
      int  len = strlen(eq->varname) + strlen("_normal_boundary_flux") + 2;

      BFT_MALLOC(bdy_flux_name, len, char);
      sprintf(bdy_flux_name, "%s_normal_boundary_flux", eq->varname);

      /* If a scalar: the scalar diffusive flux across the boundary
       *  If a vector: the vector dot the normal of the boundary face */
      int  flx_dim = (eqp->dim > 5) ? 3 : 1;
      cs_field_t  *bdy_flux_fld = cs_field_find_or_create(bdy_flux_name,
                                                          0, /* field_mask */
                                                          location_id,
                                                          flx_dim,
                                                          previous);

      eq->boundary_flux_id = cs_field_id_by_name(bdy_flux_name);

      const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
      cs_field_set_key_int(bdy_flux_fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(bdy_flux_fld, cs_field_key_id("post_vis"),
                           post_flag);

      BFT_FREE(bdy_flux_name);

    } /* Create a boundary flux field */

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the builder of the algebraic system.
 *         Set the initialize condition to all variable fields associated to
 *         each cs_equation_t structure.
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

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    const cs_equation_param_t  *eqp = eq->param;

    /* Allocate and initialize a system builder */
    eq->builder = cs_equation_init_builder(eqp, mesh);
    eq->scheme_context = eq->init_context(eqp,
                                          eq->field_id,
                                          eq->boundary_flux_id,
                                          eq->builder);

    /* Assign an initial value for the variable fields */
    if (ts->nt_cur < 1)
      eq->init_field_values(ts->t_cur, eq->field_id,
                            mesh, eqp, eq->builder, eq->scheme_context);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  }  /* Loop on equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system for this equation
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *mesh,
                         cs_equation_t              *eq)
{
  assert(eq != NULL);

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
  eq->build_system(mesh, fld->val,
                   eq->param,
                   eq->builder,
                   eq->scheme_context,
                   eq->rhs,
                   eq->matrix);

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
cs_equation_solve_deprecated(cs_equation_t   *eq)
{
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);
  cs_field_t  *fld = cs_field_by_id(eq->field_id);
  cs_real_t  *x = NULL, *b = NULL;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  const cs_equation_param_t  *eqp = eq->param;
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_sles_t  sles_param = eqp->sles_param;

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
                                                   sles_param.eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   b,
                                                   x,
                                                   0,      // aux. size
                                                   NULL);  // aux. buffers

  if (sles_param.verbosity > 0) {

    const cs_lnum_t  size = eq->n_sles_gather_elts;
    const cs_lnum_t  *row_index, *col_id;
    const cs_real_t  *d_val, *x_val;

    cs_matrix_get_msr_arrays(eq->matrix, &row_index, &col_id, &d_val, &x_val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_DBG > 1
    cs_dbg_dump_linear_system(eqp->name, size, CS_EQUATION_DBG,
                              x, b,
                              row_index, col_id, x_val, d_val);
#endif

    cs_gnum_t  nnz = row_index[size];
    cs_parall_counter(&nnz, 1);
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, n_iters, residual, nnz);

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

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Define the new field value for the current time */
  eq->update_field(x, eq->rhs, eq->param,
                   eq->builder, eq->scheme_context, fld->val);

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
 * \brief  Build and then solve the linear system for this equation when the
 *         goal is to find the steady state
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_steady_state(const cs_mesh_t            *mesh,
                               cs_equation_t              *eq)
{
  assert(eq != NULL);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
     The linear solver is called inside and the field value is updated inside
  */
  eq->solve(mesh, eq->field_id, eq->param, eq->builder, eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for an equation with an
 *         unsteady term
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(const cs_mesh_t            *mesh,
                  cs_equation_t              *eq)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty equation structure", __func__);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
     The linear solver is called inside and the field value is updated inside
  */
  eq->solve(mesh, eq->field_id, eq->param, eq->builder, eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve the related cellwise builder
 *         structures: cs_cell_builder_t and cs_cell_system_t structures
 *
 * \param[in]   eq       pointer to a \ref cs_equation_t structure
 * \param[out]  cb       pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  csys     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_get_cellwise_builders(const cs_equation_t    *eq,
                                  cs_cell_sys_t         **csys,
                                  cs_cell_builder_t     **cb)
{
  *csys = NULL;
  *cb = NULL;

  if (eq == NULL)
    return;

  if (eq->get_cw_build_structures != NULL)
    eq->get_cw_build_structures(csys, cb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         face of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 *
 * \return a pointer to an array of face values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *f_values = NULL;
  if (eq->get_face_values != NULL)
    f_values = eq->get_face_values(eq->scheme_context);

  return f_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         cell of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 *
 * \return a pointer to an array of cell values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *c_values = NULL;
  if (eq->get_cell_values != NULL)
    c_values = eq->get_cell_values(eq->scheme_context);

  return c_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         vertex of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 *
 * \return a pointer to an array of vertex values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_vertex_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *v_values = NULL;
  if (eq->get_vertex_values != NULL)
    v_values = eq->get_vertex_values(eq->scheme_context);

  return v_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive flux across all boundary faces
 *         According to the space discretization scheme, the size of the
 *         resulting array differs.
 *         For Vb and VCb schemes, this array relies on the bf2v adjacency.
 *
 * \param[in]      t_eval     time at which one performs the property evaluation
 * \param[in]      eq         pointer to a cs_equation_t structure
 * \param[in, out] diff_flux  value of the diffusive part of the flux
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_boundary_diff_flux(cs_real_t              t_eval,
                                       const cs_equation_t   *eq,
                                       cs_real_t             *diff_flux)
{
  if (diff_flux == NULL)
    return;

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq, __func__);

  const cs_equation_param_t  *eqp = eq->param;

  assert(eqp != NULL);
  if (eqp->dim > 1)
    bft_error(__FILE__, __LINE__, 0, "%s: (Eq. %s) Not implemented",
              __func__, eqp->name);

  /* Scalar-valued equation */
  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *p_v = cs_equation_get_vertex_values(eq);

      cs_cdovb_scaleq_boundary_diff_flux(t_eval,
                                         eqp,
                                         p_v,
                                         eq->builder,
                                         diff_flux);
    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *p_v = cs_equation_get_vertex_values(eq);
      const cs_real_t  *p_c = cs_equation_get_cell_values(eq);

      cs_cdovcb_scaleq_boundary_diff_flux(t_eval,
                                          eqp,
                                          p_v,
                                          p_c,
                                          eq->builder,
                                          diff_flux);
    }
    break;

  case CS_SPACE_SCHEME_CDOFB:
    {
      const cs_real_t  *p_f = cs_equation_get_face_values(eq);
      const cs_real_t  *p_c = cs_equation_get_cell_values(eq);

      cs_cdofb_scaleq_boundary_diff_flux(t_eval,
                                         eqp,
                                         p_f,
                                         p_c,
                                         eq->builder,
                                         diff_flux);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: (Eq. %s). Not implemented.", __func__, eqp->name);

  } /* End of switch */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a plane defined
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
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq, __func__);

  /* Get the mesh location id from its name */
  const int  ml_id = cs_mesh_location_get_id_by_name(ml_name);
  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid mesh location name %s.\n"
              " This mesh location is not already defined.\n",
              __func__, ml_name);

  const char  emsg[] = " %s: Computation of the diffusive and convective flux"
    " across a plane\n is not available for equation %s\n";

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  cs_equation_param_t  *eqp = eq->param;
  assert(eqp != NULL && fld != NULL);

  if (eqp->dim > 1)
    bft_error(__FILE__, __LINE__, 0, emsg, __func__, eqp->name);


  /* Scalar-valued equation */
  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    cs_cdovb_scaleq_flux_across_plane(direction,
                                      fld->val,
                                      eqp,
                                      ml_id,
                                      eq->builder,
                                      eq->scheme_context,
                                      diff_flux, conv_flux);
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    cs_cdovcb_scaleq_flux_across_plane(direction,
                                       fld->val,
                                       eqp,
                                       ml_id,
                                       eq->builder,
                                       eq->scheme_context,
                                       diff_flux, conv_flux);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, emsg, __func__, eqp->name);

  } /* End of switch */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all cell faces.
 *         Primal or dual faces are considered according to the space scheme.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      location    indicate where the flux has to be computed
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diff_flux_cellwise(const cs_equation_t   *eq,
                                       cs_flag_t              location,
                                       cs_real_t              t_eval,
                                       cs_real_t             *diff_flux)
{
  if (diff_flux == NULL)
    return;
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq, __func__);

  cs_equation_param_t  *eqp = eq->param;
  assert(eqp != NULL);

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  const char  fmsg[] = " %s: (Eq. %s) Stop computing the diffusive flux.\n"
    " This functionality is not available for this scheme.";
  const char  lmsg[] = " %s: (Eq. %s) Stop computing the diffusive flux.\n"
    " This mesh location is not available for this scheme.";

  if (eqp->dim > 1)
    bft_error(__FILE__, __LINE__, 0, fmsg, __func__, eqp->name);

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    if (cs_flag_test(location, cs_flag_primal_cell))
      cs_cdovb_scaleq_diff_flux_in_cells(fld->val,
                                         eqp,
                                         t_eval,
                                         eq->builder,
                                         eq->scheme_context,
                                         diff_flux);
    else if (cs_flag_test(location, cs_flag_dual_face_byc))
      cs_cdovb_scaleq_diff_flux_dfaces(fld->val,
                                       eqp,
                                       t_eval,
                                       eq->builder,
                                       eq->scheme_context,
                                       diff_flux);
    else
      bft_error(__FILE__, __LINE__, 0, lmsg, __func__, eqp->name);
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    if (cs_flag_test(location, cs_flag_primal_cell))
      cs_cdovcb_scaleq_diff_flux_in_cells(fld->val,
                                          eqp,
                                          t_eval,
                                          eq->builder,
                                          eq->scheme_context,
                                          diff_flux);
    else if (cs_flag_test(location, cs_flag_dual_face_byc))
      cs_cdovcb_scaleq_diff_flux_dfaces(fld->val,
                                        eqp,
                                        t_eval,
                                        eq->builder,
                                        eq->scheme_context,
                                        diff_flux);
    else
      bft_error(__FILE__, __LINE__, 0, lmsg, __func__, eqp->name);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, fmsg, __func__, eqp->name);

  } /* End of switch on the space scheme */

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
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq, __func__);
  assert(v_gradient != NULL);

  const cs_equation_param_t  *eqp = eq->param;

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVCB:
    cs_cdovcb_scaleq_vtx_gradient(fld->val, eq->builder, eq->scheme_context,
                                  v_gradient);
    break;

  case CS_SPACE_SCHEME_CDOVB:   /* --> wall distance */
  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of scheme for equation %s when computing"
              " the gradient at vertices", __func__, eqp->name);
    break;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute and post-process Peclet number if requested
 *
 * \param[in]      eq       pointer to a cs_equation_t structure
 * \param[in]      ts       pointer to a cs_time_step_t struct.
 * \param[in, out] peclet   pointer to an array storing the resulting Peclet
 *                          number in each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_peclet(const cs_equation_t        *eq,
                           const cs_time_step_t       *ts,
                           cs_real_t                   peclet[])
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq, __func__);
  assert(peclet != NULL);

  const cs_equation_param_t  *eqp = eq->param;

  /* Check if the computation of the Peclet number is requested */
  if (!(eqp->process_flag & CS_EQUATION_POST_PECLET))
    return;

  if (eqp->diffusion_property == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Computation of the Peclet number is requested for\n"
              " equation %s but no diffusion property is set.\n",
              __func__, eqp->name);
  if (eqp->adv_field == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Computation of the Peclet number is requested for\n"
              " equation %s but no advection field is set.\n",
              __func__, eqp->name);

  if (eq->main_ts_id > -1)    /* Activate timer statistics */
    cs_timer_stats_start(eq->main_ts_id);

  /* Compute the Peclet number in each cell */
  cs_advection_get_peclet(eqp->adv_field,
                          eqp->diffusion_property,
                          ts->t_cur,
                          peclet);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write into the restart file additionnal arrays (not defined as
 *         fields) but useful for the checkpoint/restart process
 *
 * \param[in, out]  restart    pointer to a \ref cs_restart_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_read_extra_restart(cs_restart_t   *restart)
{
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */

    if (eq->read_restart != NULL)
      eq->read_restart(restart, eq->param->name, eq->scheme_context);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write into the restart file additionnal arrays (not defined as
 *         fields) but useful for the checkpoint/restart process
 *
 * \param[in, out]  restart    pointer to a \ref cs_restart_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_write_extra_restart(cs_restart_t   *restart)
{
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */

    if (eq->write_restart != NULL)
      eq->write_restart(restart, eq->param->name, eq->scheme_context);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to all equations
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_post_balance(const cs_mesh_t            *mesh,
                         const cs_cdo_connect_t     *connect,
                         const cs_cdo_quantities_t  *cdoq,
                         const cs_time_step_t       *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(cdoq);

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */
    const cs_equation_param_t  *eqp = eq->param;

    /* Check if the computation of the balance is requested */
    if (!(eqp->process_flag & CS_EQUATION_POST_BALANCE))
      continue;

    if (eq->compute_balance != NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Balance for equation %s is requested but\n"
                " this functionality is not available yet.\n",
                __func__, eqp->name);

    if (eq->main_ts_id > -1)    /* Activate timer statistics */
      cs_timer_stats_start(eq->main_ts_id);

    cs_equation_balance_t  *b = eq->compute_balance(eqp,
                                                    eq->builder,
                                                    eq->scheme_context);

    char *postlabel = NULL;
    int len = strlen(eqp->name) + 13 + 1;
    BFT_MALLOC(postlabel, len, char);

    switch (eqp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      {
        sprintf(postlabel, "%s.Balance", eqp->name);

        cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                                 CS_POST_WRITER_DEFAULT,
                                 postlabel,
                                 eqp->dim,
                                 false,
                                 false,
                                 CS_POST_TYPE_cs_real_t,
                                 b->balance,
                                 ts);

        if (cs_equation_param_has_diffusion(eqp))
          _post_balance_at_vertices(eq, ts, "Diff", postlabel,
                                    b->diffusion_term);

        if (cs_equation_param_has_convection(eqp))
          _post_balance_at_vertices(eq, ts, "Adv", postlabel,
                                    b->advection_term);

        if (cs_equation_param_has_time(eqp))
          _post_balance_at_vertices(eq, ts, "Time", postlabel,
                                    b->unsteady_term);

        if (cs_equation_param_has_reaction(eqp))
          _post_balance_at_vertices(eq, ts, "Reac", postlabel,
                                    b->reaction_term);

        if (cs_equation_param_has_sourceterm(eqp))
          _post_balance_at_vertices(eq, ts, "Src", postlabel,
                                    b->source_term);

      }
      break;

    default:
      break;
    }

    sprintf(postlabel, "%s.BdyFlux", eqp->name);

    /* Post-process the boundary fluxes (diffusive and convective) */
    cs_post_write_var(CS_POST_MESH_BOUNDARY,
                      CS_POST_WRITER_DEFAULT,
                      postlabel,
                      1,
                      true,             // interlace
                      true,             // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      NULL,             // values on cells
                      NULL,             // values at internal faces
                      b->boundary_term, // values at border faces
                      ts);              // time step management struct.

    /* Free buffers */
    BFT_FREE(postlabel);
    cs_equation_balance_destroy(&b);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to equations according to the
 *         type of numerical scheme (for the space discretization)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_post(void)
{
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */
    const cs_equation_param_t  *eqp = eq->param;

    if (eq->main_ts_id > -1)    /* Activate timer statistics */
      cs_timer_stats_start(eq->main_ts_id);
    assert(eq->postprocess != NULL);

    /* Perform post-processing specific to a numerical scheme */
    eq->postprocess(eqp->name,
                    cs_field_by_id(eq->field_id),
                    eqp,
                    eq->builder,
                    eq->scheme_context);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
