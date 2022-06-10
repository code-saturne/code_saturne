/*============================================================================
 * Functions to handle the cs_equation_t structure and its related structures
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_cdo_toolbox.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovb_vecteq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdoeb_vecteq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_equation_bc.h"
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
 * Local macro definitions
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
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
solve_steady_state_cdo_equation(const char       *eqname);

void
solve_cdo_equation(bool         cur2prev,
                   const char  *eqname);

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
 *         Deprecated function.
 *
 * \param[in, out] eq_to_cast   pointer to generic builder structure
 * \param[in, out] p_x          pointer of pointer to the solution array
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_fb_solving(void              *eq_to_cast,
                    cs_real_t         *p_x[])
{
  cs_equation_t  *eq = (cs_equation_t  *)eq_to_cast;
  cs_equation_param_t  *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;

  const cs_range_set_t  *rset = cs_equation_builder_get_range_set(eqb, 0);
  const cs_matrix_t *matrix = cs_equation_builder_get_matrix(eqb, 0);
  const cs_lnum_t  n_dofs = eqp->dim * rset->n_elts[1];
  const cs_real_t  *f_values = eq->get_face_values(eq->scheme_context, false);
  const int  stride = 1;  /* Since the global numbering is adapted in each
                             case (scalar-, vector-valued equations) */

  assert(f_values != NULL);

  cs_real_t  *x = NULL;
  BFT_MALLOC(x, CS_MAX(n_dofs, cs_matrix_get_n_columns(matrix)), cs_real_t);

  /* x and the right-hand side are a "gathered" view of field->val and the
   * right-hand side respectively through the range set operation.
   *
   *  Their size is equal to n_sles_gather_elts <= n_dofs
   */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_cdo_system_helper_t  *sh = eqb->system_helper;

    /* Compact numbering to fit the algebraic decomposition */

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,  /* type */
                        stride,        /* stride */
                        f_values,      /* in: size = n_dofs */
                        x);            /* out: size = n_sles_gather_elts */

    cs_interface_set_sum(rset->ifs,
                         n_dofs, stride, false, CS_REAL_TYPE,
                         sh->rhs);

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,  /* type */
                        stride,        /* stride */
                        sh->rhs,       /* in: size = n_dofs */
                        sh->rhs);      /* out: size = n_sles_gather_elts */

  }
  else { /* Serial mode *** without periodicity *** */

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (n_dofs > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < n_dofs; i++)
      x[i] = f_values[i];
#else
    memcpy(x, f_values, n_dofs * sizeof(cs_real_t));
#endif

  }

  /* Return pointers */

  *p_x = x;
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

  eq->postprocess = cs_hho_scaleq_extra_post;
  eq->read_restart = cs_hho_scaleq_read_restart;
  eq->write_restart = cs_hho_scaleq_write_restart;

  /* Function pointers to retrieve values at mesh locations */

  eq->get_vertex_values = NULL;
  eq->get_edge_values = NULL;
  eq->get_face_values = cs_hho_scaleq_get_face_values;
  eq->get_cell_values = cs_hho_scaleq_get_cell_values;

  /* Deprecated functions */

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

  eq->postprocess = cs_hho_vecteq_extra_post;
  eq->read_restart = cs_hho_vecteq_read_restart;
  eq->write_restart = cs_hho_vecteq_write_restart;

  /* Function pointers to retrieve values at mesh locations */

  eq->get_vertex_values = NULL;
  eq->get_edge_values = NULL;
  eq->get_face_values = cs_hho_vecteq_get_face_values;
  eq->get_cell_values = cs_hho_vecteq_get_cell_values;

  /* Deprecated functions */

  eq->build_system = cs_hho_vecteq_build_system;
  eq->prepare_solving = _prepare_fb_solving;
  eq->update_field = cs_hho_vecteq_update_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a field (type=variable) to an equation
 *
 * \param[in]       n_previous   number of previous states to keep
 * \param[in, out]  eq           pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_field(int               n_previous,
           cs_equation_t    *eq)
{
  if (eq == NULL)
    return;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  /* Redundant definition to handle C/FORTRAN */

  bool  previous;
  if (n_previous == 0)
    previous = false;
  else if (n_previous > 0)
    previous = true;
  else {

    previous = false; /* avoid a compiler warning */
    bft_error(__FILE__, __LINE__, 0,
              "%s: Expected value for n_previous is > -1. Here %d\n"
              "%s: Eq. \"%s\"\n",
              __func__, n_previous, __func__, eqp->name);

  }

  /* Associate a predefined mesh_location_id to this field */

  int  location_id = -1; /* initialize values to avoid a warning */

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
              "%s: Space scheme for eq. \"%s\" is incompatible with a field.\n"
              "%s: Stop adding a cs_field_t structure.\n",
              __func__, eqp->name, __func__);
    break;
  }

  if (location_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid mesh location id (= -1) for the field associated"
              " to Eq. \"%s\"\n", __func__, eqp->name);

  /* Store the related field id */

  eq->field_id = cs_variable_cdo_field_create(eq->varname,
                                              NULL, /* label */
                                              location_id,
                                              eqp->dim,
                                              n_previous);

  /* SLES is associated to a field_id */

  eqp->sles_param->field_id = eq->field_id;

  if (eqp->post_flag & CS_EQUATION_POST_NORMAL_FLUX) {

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
    assert(eqp != NULL && eqp->name != NULL);
    if (strlen(eqp->name) == len_in)
      if (strcmp(eqname, eqp->name) == 0)
        return _eq;

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the pointer to a cs_equation_t structure thanks to the field
 *         name of the variable field associated to a cs_equation_t structure
 *
 * \param[in]  field_name       name of the field
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_field_name(const char    *field_name)
{
  if (field_name == NULL)
    return NULL;

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    assert(eq != NULL);

    if (cs_equation_has_field_name(eq, field_name))
      return eq;

  } /* Loop on equations */

  return NULL;
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
 *         cs_equation_t structure based on the equation name
 *
 * If no equation matches the given name but a field does, equation
 * parameter structure associated to the field will be returned instead.
 * This allows using this function with non-CDO (legacy) fields.
 *
 * \param[in]  eqname       name of the equation
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_by_name(const char    *eqname)
{
  cs_equation_param_t *eq_param = NULL;

  if (eqname != NULL) {

    cs_equation_t  *eq = cs_equation_by_name(eqname);
    if (eq != NULL)
      eq_param = eq->param;

    else {
      cs_field_t *f = cs_field_by_name_try(eqname);
      if (f != NULL)
        eq_param = cs_field_get_equation_param(f);
    }

  }

  return eq_param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure related to a
 *         cs_equation_t structure thanks to the field name of the variable
 *         field associated to a cs_equation_t structure
 *
 * \param[in]  field_name       name of the field
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_by_field_name(const char    *field_name)
{
  if (field_name == NULL)
    return NULL;

  cs_equation_t  *eq = cs_equation_by_field_name(field_name);

  if (eq == NULL)
    return NULL;
  else
    return eq->param;
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
 * \brief  Return the range set structure associated to a cs_equation_t
 *         structure. One assumes that there is only one block (it could be a
 *         split block) otherwise this means that one handles systems of
 *         equations.
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_range_set_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const cs_range_set_t *
cs_equation_get_range_set(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  else { /* Equation has been allocated */

    cs_equation_builder_t  *eqb = eq->builder;

    if (eqb == NULL)
      return NULL;
    else {

      return cs_equation_builder_get_range_set(eqb, 0);
    }

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the global number of degrees of freedom associated to this
 *         cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return a global number of degrees of freedom (DoFs)
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_equation_get_global_n_dofs(const cs_equation_t         *eq,
                              const cs_cdo_quantities_t   *cdoq)
{
  if (eq == NULL || cdoq == NULL)
    return 0;

  switch (cs_equation_get_space_scheme(eq)) {

  case CS_SPACE_SCHEME_CDOVB:
    if (cs_glob_n_ranks > 1)
      return cdoq->n_g_vertices;
    else
      return (cs_gnum_t)cdoq->n_vertices;
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    if (cs_glob_n_ranks > 1)
      return cdoq->n_g_vertices + cdoq->n_g_cells;
    else
      return (cs_gnum_t)(cdoq->n_vertices + cdoq->n_cells);
    break;

  case CS_SPACE_SCHEME_CDOEB:
    if (cs_glob_n_ranks > 1)
      return cdoq->n_g_edges;
    else
      return (cs_gnum_t)(cdoq->n_edges);
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    if (cs_glob_n_ranks > 1)
      return cdoq->n_g_faces + cdoq->n_g_cells;
    else
      return (cs_gnum_t)(cdoq->n_faces + cdoq->n_cells);
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    if (cs_glob_n_ranks > 1)
      return CS_N_DOFS_FACE_1ST*cdoq->n_g_faces
        + CS_N_DOFS_CELL_1ST*cdoq->n_g_cells;
    else
      return (cs_gnum_t)(CS_N_DOFS_FACE_1ST*cdoq->n_faces
                         + CS_N_DOFS_CELL_1ST*cdoq->n_cells);
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    if (cs_glob_n_ranks > 1)
      return CS_N_DOFS_FACE_2ND*cdoq->n_g_faces
        + CS_N_DOFS_CELL_2ND*cdoq->n_g_cells;
    else
      return (cs_gnum_t)(CS_N_DOFS_FACE_2ND*cdoq->n_faces
                         + CS_N_DOFS_CELL_2ND*cdoq->n_cells);
    break;

  default:
    return 0;
  }

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
 * \brief  Add a hook function to enable an advanced control during the
 *         cellwise system building.
 *         Only for an advanced usage. The context may be set to NULL if there
 *         is no need to get additional information.
 *
 * \param[in, out] eq        pointer to the cs_equation_t stucture to update
 * \param[in]      context   pointer to a structure for additional information
 * \param[in]      func      pointer to the user function
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_build_hook(cs_equation_t               *eq,
                           void                        *context,
                           cs_equation_build_hook_t    *func)
{
  if (eq == NULL)
    return;

  cs_equation_param_t  *eqp = eq->param;
  assert(eqp != NULL);

  if (eq->builder == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Initialization of equation %s has not been done yet.\n"
              " Please call this operation later in"
              " cs_user_extra_operations_initialize() for instance.",
              __func__, eqp->name);

  cs_equation_builder_t   *eqb = eq->builder;

  eqb->hook_context = context;
  eqb->hook_function = func;
  eqp->flag |= CS_EQUATION_BUILD_HOOK;

  /* Add an entry in the setup log file (this is done after the main setup
   * log but one needs to initialize equations before calling this function) */

  cs_log_printf(CS_LOG_SETUP, " Equation %s: Add a user hook function\n",
                eqp->name);
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
 * \brief  Return a pointer to a structure useful to handle low-level
 *         operations for the given equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  structure storing the main structure associated to an equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_core_t
cs_equation_get_core(const cs_equation_t    *eq)
{
  cs_equation_core_t  core =
    { .param = NULL, .builder = NULL, .scheme_context = NULL };

  if (eq != NULL) {

    core.param = eq->param;
    core.builder = eq->builder;
    core.scheme_context = eq->scheme_context;

  }

  return core;
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
    case CS_TIME_SCHEME_BDF2:
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
        (eq->param->space_scheme == CS_SPACE_SCHEME_CDOFB) ||
        (eq->param->space_scheme == CS_SPACE_SCHEME_CDOEB))
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

  case CS_EQUATION_TYPE_GROUNDWATER:
  case CS_EQUATION_TYPE_MAXWELL:
  case CS_EQUATION_TYPE_NAVSTO:
  case CS_EQUATION_TYPE_PREDEFINED:
  case CS_EQUATION_TYPE_SOLIDIFICATION:
  case CS_EQUATION_TYPE_THERMAL:
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

  size_t  len = strlen(varname);
  BFT_MALLOC(eq->varname, len + 1, char);
  strncpy(eq->varname, varname, len);
  eq->varname[len] = '\0';

  eq->param = cs_equation_param_create(eqname, eqtype, dim, default_bc);

  eq->field_id = -1;           /* This field is created in a second step */
  eq->boundary_flux_id = -1;   /* Not always defined (done in a second step) */

  /* Builder structure for this equation */

  eq->builder = NULL;
  eq->scheme_context = NULL;

  /* Pointers of function */

  eq->init_context = NULL;
  eq->free_context = NULL;

  /* Extra-operations */

  eq->compute_balance = NULL;
  eq->apply_stiffness = NULL;
  eq->postprocess = NULL;
  eq->current_to_previous = NULL;

  /* Restart */

  eq->read_restart = NULL;
  eq->write_restart = NULL;

  /* Function pointers to retrieve values at mesh locations */

  eq->get_vertex_values = NULL;
  eq->get_edge_values = NULL;
  eq->get_face_values = NULL;
  eq->get_cell_values = NULL;

  /* New functions */

  eq->init_field_values = NULL;
  eq->solve = NULL;
  eq->solve_steady_state = NULL;

  /* Deprecated functions */

  eq->set_dir_bc = NULL;
  eq->build_system = NULL;
  eq->update_field = NULL;

  /* Set timer statistic structure to a default value */

  eq->main_ts_id = cs_timer_stats_id_by_name(eqname);
  if (eq->main_ts_id < 0)
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
 * \brief  Add a new user transport equation and set a first set of parameters
 *         If time_pty is NULL, then no unsteady term is added.
 *         If adv is NULL, then no advection term is added.
 *         If diff_pty is NULL, then no diffusion term is added.
 *
 * \param[in] eqname       name of the equation
 * \param[in] varname      name of the variable associated to this equation
 * \param[in] dim          dimension of the unknow attached to this equation
 * \param[in] default_bc   type of boundary condition set by default
 * \param[in] time_pty     property related to the unsteady term
 * \param[in] adv          advection field
 * \param[in] diff_pty     property related to the diffusion term
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add_user_tracer(const char            *eqname,
                            const char            *varname,
                            int                    dim,
                            cs_param_bc_type_t     default_bc,
                            cs_property_t         *time_pty,
                            cs_adv_field_t        *adv,
                            cs_property_t         *diff_pty)
{
  cs_equation_t  *eq = cs_equation_add_user(eqname, varname, dim, default_bc);
  assert(eq != NULL);

  /* Add an advection term */

  if (adv != NULL)
    cs_equation_add_advection(eq->param, adv);

  /* Add the unsteady term */

  if (time_pty != NULL)
    cs_equation_add_time(eq->param, time_pty);

  /* Add the diffusion term */

  if (diff_pty != NULL)
    cs_equation_add_diffusion(eq->param, diff_pty);

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

    /* Free the associated builder structure */

    cs_equation_builder_free(&(eq->builder));
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
  /* Check if there is something to output */

  if (_n_equations < 1)
    return;

  bool  output = false;
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    const cs_equation_param_t  *eqp = eq->param;

    if (eqp->flag & CS_EQUATION_INSIDE_SYSTEM)
      continue;
    else
      output = true;

  }

  if (!output)
    return;

  cs_log_printf(CS_LOG_PERFORMANCE, "\n%-36s %9s %9s %9s\n",
                " ", "Build", "Solve", "Extra");

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];

    /* Display high-level timer counter related to the current equation before
       deleting the structure */

    cs_equation_builder_log_performance(eq->param, eq->builder);

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
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  for (int  eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    assert(eq->param != NULL);

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    cs_log_printf(CS_LOG_SETUP,
                  "Summary of settings for \"%s\" eq. (variable: \"%s\")\n",
                  eq->param->name, eq->varname);
    cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h2);

    cs_equation_param_log(eq->param);

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

    cs_equation_param_set(eq->param, key, keyval);

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
 * \brief  Set shared pointers to the main structures. Associate these
 *         structures among the activated class of discretization schemes.
 *         Initialize
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
cs_equation_init_sharing(const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *time_step,
                         cs_flag_t                    eb_scheme_flag,
                         cs_flag_t                    fb_scheme_flag,
                         cs_flag_t                    vb_scheme_flag,
                         cs_flag_t                    vcb_scheme_flag,
                         cs_flag_t                    hho_scheme_flag)
{
  if (vb_scheme_flag > 0 || vcb_scheme_flag > 0) {

    if (vb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
      cs_cdovb_scaleq_init_sharing(quant, connect, time_step);

    if (vcb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
      cs_cdovcb_scaleq_init_sharing(quant, connect, time_step);

    if (vb_scheme_flag & CS_FLAG_SCHEME_VECTOR)
      cs_cdovb_vecteq_init_sharing(quant, connect, time_step);

  } /* Vertex-based class of discretization schemes */

  if (eb_scheme_flag > 0) {

      /* This is a vector-valued equation but the DoF is scalar-valued since
       * it is a circulation associated to each edge */

    if (eb_scheme_flag  & CS_FLAG_SCHEME_SCALAR)
      cs_cdoeb_vecteq_init_sharing(quant, connect, time_step);

  } /* Edge-based class of discretization schemes */

  if (fb_scheme_flag > 0) {

    if (cs_flag_test(fb_scheme_flag,
                     CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR))
      cs_cdofb_scaleq_init_sharing(quant, connect, time_step);

    if (cs_flag_test(fb_scheme_flag,
                     CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR))
      cs_cdofb_vecteq_init_sharing(quant, connect, time_step);

  } /* Face-based class of discretization schemes */

  if (hho_scheme_flag > 0) {

    if (hho_scheme_flag & CS_FLAG_SCHEME_SCALAR)
      cs_hho_scaleq_init_sharing(hho_scheme_flag, quant, connect, time_step);

    if (hho_scheme_flag & CS_FLAG_SCHEME_VECTOR)
      cs_hho_vecteq_init_sharing(hho_scheme_flag, quant, connect, time_step);

  } /* Higher-order face-based class of discretization schemes */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared local structures among the discretization schemes
 *
 * \param[in]  vb_scheme_flag   metadata for Vb schemes
 * \param[in]  vcb_scheme_flag  metadata for V+C schemes
 * \param[in]  eb_scheme_flag   metadata for Eb schemes
 * \param[in]  fb_scheme_flag   metadata for Fb schemes
 * \param[in]  hho_scheme_flag  metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_finalize_sharing(cs_flag_t    vb_scheme_flag,
                             cs_flag_t    vcb_scheme_flag,
                             cs_flag_t    eb_scheme_flag,
                             cs_flag_t    fb_scheme_flag,
                             cs_flag_t    hho_scheme_flag)
{
  /* Free common local structures specific to a numerical scheme */

  if (vb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdovb_scaleq_finalize_sharing();

  if (vb_scheme_flag & CS_FLAG_SCHEME_VECTOR)
    cs_cdovb_vecteq_finalize_sharing();

  if (vcb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdovcb_scaleq_finalize_sharing();

  if (eb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdoeb_vecteq_finalize_sharing();

  if (fb_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_cdofb_scaleq_finalize_sharing();

  if (fb_scheme_flag & CS_FLAG_SCHEME_VECTOR)
    cs_cdofb_vecteq_finalize_sharing();

  if (hho_scheme_flag & CS_FLAG_SCHEME_SCALAR)
    cs_hho_scaleq_finalize_sharing();

  if (hho_scheme_flag & CS_FLAG_SCHEME_VECTOR)
    cs_hho_vecteq_finalize_sharing();
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
      cs_equation_param_set(eqp, CS_EQKEY_TIME_SCHEME, "steady");

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
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */

        if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE) {

          eq->solve_steady_state = cs_cdovb_scaleq_solve_steady_state_incr;

          switch (eqp->time_scheme) {
          case CS_TIME_SCHEME_STEADY:
            eq->solve = eq->solve_steady_state;
            break;
          case CS_TIME_SCHEME_EULER_IMPLICIT:
            eq->solve = cs_cdovb_scaleq_solve_implicit_incr;
            break;
          case CS_TIME_SCHEME_THETA:
          case CS_TIME_SCHEME_CRANKNICO:
          case CS_TIME_SCHEME_BDF2:
          default:
            bft_error(__FILE__, __LINE__, 0,
                      "%s: Eq. %s. This time scheme is not yet implemented",
                      __func__, eqp->name);
          }

        }
        else {

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

          case CS_TIME_SCHEME_BDF2:
          default:
            bft_error(__FILE__, __LINE__, 0,
                      "%s: Eq. %s. This time scheme is not yet implemented",
                      __func__, eqp->name);
          }

        } /* Incremental solve or not */

        eq->compute_balance = cs_cdovb_scaleq_balance;
        eq->apply_stiffness = cs_cdovb_scaleq_apply_stiffness;
        eq->postprocess = cs_cdovb_scaleq_extra_post;
        eq->current_to_previous = cs_cdovb_scaleq_current_to_previous;

        eq->read_restart = NULL;
        eq->write_restart = NULL;

        /* Function pointers to retrieve values at mesh locations */

        eq->get_vertex_values = cs_cdovb_scaleq_get_vertex_values;
        eq->get_edge_values = NULL;
        eq->get_face_values = NULL;
        eq->get_cell_values = cs_cdovb_scaleq_get_cell_values;

        eq->get_cw_build_structures = cs_cdovb_scaleq_get;

      }
      else if (eqp->dim == 3) {

        eq->init_context = cs_cdovb_vecteq_init_context;
        eq->free_context = cs_cdovb_vecteq_free_context;
        eq->init_field_values = cs_cdovb_vecteq_init_values;

        /* Deprecated */

        eq->set_dir_bc = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */

        if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. Incremental form is not available.\n",
                    __func__, eqp->name);

        eq->solve_steady_state = cs_cdovb_vecteq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
        eq->solve = eq->solve_steady_state;
          break;
        case CS_TIME_SCHEME_BDF2:
        case CS_TIME_SCHEME_EULER_IMPLICIT:
        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdovb_vecteq_extra_post;
        eq->current_to_previous = cs_cdovb_vecteq_current_to_previous;

        eq->read_restart = NULL;
        eq->write_restart = NULL;

        /* Function pointers to retrieve values at mesh locations */

        eq->get_vertex_values = cs_cdovb_vecteq_get_vertex_values;
        eq->get_edge_values = NULL;
        eq->get_face_values = NULL;
        eq->get_cell_values = cs_cdovb_vecteq_get_cell_values;

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
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */

        if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. Incremental form is not available.\n",
                    __func__, eqp->name);

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

        case CS_TIME_SCHEME_BDF2:
        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdovcb_scaleq_extra_post;
        eq->current_to_previous = cs_cdovcb_scaleq_current_to_previous;

        eq->read_restart = cs_cdovcb_scaleq_read_restart;
        eq->write_restart = cs_cdovcb_scaleq_write_restart;

        /* Function pointers to retrieve values at mesh locations */

        eq->get_vertex_values = cs_cdovcb_scaleq_get_vertex_values;
        eq->get_edge_values = NULL;
        eq->get_face_values = NULL;
        eq->get_cell_values = cs_cdovcb_scaleq_get_cell_values;

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
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */

        if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. Incremental form is not available.\n",
                    __func__, eqp->name);

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

        case CS_TIME_SCHEME_BDF2:
        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->compute_balance = cs_cdofb_scaleq_balance;
        eq->apply_stiffness = NULL;
        eq->postprocess = cs_cdofb_scaleq_extra_post;
        eq->current_to_previous = cs_cdofb_scaleq_current_to_previous;

        eq->read_restart = cs_cdofb_scaleq_read_restart;
        eq->write_restart = cs_cdofb_scaleq_write_restart;

        /* Function pointers to retrieve values at mesh locations */

        eq->get_vertex_values = NULL;
        eq->get_edge_values = NULL;
        eq->get_face_values = cs_cdofb_scaleq_get_face_values;
        eq->get_cell_values = cs_cdofb_scaleq_get_cell_values;

        eq->get_cw_build_structures = cs_cdofb_scaleq_get;

      }
      else if (eqp->dim == 3) {

        eq->init_context = cs_cdofb_vecteq_init_context;
        eq->free_context = cs_cdofb_vecteq_free_context;
        eq->init_field_values = cs_cdofb_vecteq_init_values;

        /* Deprecated */

        eq->set_dir_bc = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */

        if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. Incremental form is not available.\n",
                    __func__, eqp->name);

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

        case CS_TIME_SCHEME_BDF2:
          eq->solve = NULL; /* cs_cdofb_vecteq_solve_bdf2 */
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdofb_vecteq_extra_post;
        eq->current_to_previous = cs_cdofb_vecteq_current_to_previous;

        eq->read_restart = cs_cdofb_vecteq_read_restart;
        eq->write_restart = cs_cdofb_vecteq_write_restart;

        /* Function pointers to retrieve values at mesh locations */

        eq->get_vertex_values = NULL;
        eq->get_edge_values = NULL;
        eq->get_face_values = cs_cdofb_vecteq_get_face_values;
        eq->get_cell_values = cs_cdofb_vecteq_get_cell_values;

        eq->get_cw_build_structures = cs_cdofb_vecteq_get;
      }
      else
        bft_error(__FILE__, __LINE__, 0, sv_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_CDOEB:
      if (eqp->dim == 3) {

        /* This is a vector-valued equation but the DoF is scalar-valued since
         * it is a circulation associated to each edge */

        eq->init_context = cs_cdoeb_vecteq_init_context;
        eq->free_context = cs_cdoeb_vecteq_free_context;
        eq->init_field_values = cs_cdoeb_vecteq_init_values;

        /* Deprecated */

        eq->set_dir_bc = NULL;
        eq->build_system = NULL;
        eq->prepare_solving = NULL;
        eq->update_field = NULL;

        /* New mechanism */

        if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. Incremental form is not available.\n",
                    __func__, eqp->name);

        eq->solve_steady_state = cs_cdoeb_vecteq_solve_steady_state;
        switch (eqp->time_scheme) {
        case CS_TIME_SCHEME_STEADY:
          eq->solve = eq->solve_steady_state;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Eq. %s. This time scheme is not yet implemented",
                    __func__, eqp->name);
        }

        eq->postprocess = cs_cdoeb_vecteq_extra_post;
        eq->current_to_previous = cs_cdoeb_vecteq_current_to_previous;

        eq->read_restart = cs_cdoeb_vecteq_read_restart;
        eq->write_restart = cs_cdoeb_vecteq_write_restart;

        /* Function pointers to retrieve values at mesh locations */

        eq->get_vertex_values = NULL;
        eq->get_edge_values = cs_cdoeb_vecteq_get_edge_values;
        eq->get_face_values = NULL;
        eq->get_cell_values = cs_cdoeb_vecteq_get_cell_values;

        eq->get_cw_build_structures = cs_cdoeb_vecteq_get;

      }
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P0:
      if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Eq. %s. Incremental form is not available.\n",
                  __func__, eqp->name);

      if (eqp->dim == 1) /* Set pointers of function */
        _set_scal_hho_function_pointers(eq);
      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P1:
      if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Eq. %s. Incremental form is not available.\n",
                  __func__, eqp->name);

      if (eqp->dim == 1) /* Set pointers of function */
        _set_scal_hho_function_pointers(eq);

      else if (eqp->dim == 3) /* Set pointers of function */
        _set_vect_hho_function_pointers(eq);

      else
        bft_error(__FILE__, __LINE__, 0, s_err_msg, __func__);

      break;

    case CS_SPACE_SCHEME_HHO_P2:
      if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Eq. %s. Incremental form is not available.\n",
                  __func__, eqp->name);

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
 * \brief  Create a field structure related to the predefined equation given as
 *         parameter. This includes an equation associated to all modules and
 *         also wall distance or mesh deformation for instance
 *
 *         When an automatic behavior is asked then one checks the flag
 *         CS_EQUATION_UNSTEADY to decide. One can force the behavior when
 *         handling predefined equations since more complex situations can
 *         occur such as a steady computation with non-linearities (in which
 *         case one wants a field with a previous state)
 *
 * \param[in]       n_previous     number of previous states to keep
 *                                 -1 means automatic
 * \param[in, out]  eq             pointer to an equation structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_predefined_create_field(int               n_previous,
                                    cs_equation_t    *eq)
{
  if (eq == NULL)
    return;

  cs_equation_param_t  *eqp = eq->param;

  if (eqp->type == CS_EQUATION_TYPE_USER)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Only predefined equation are managed with this function.\n"
              "%s: Eq. \"%s\"\n",
              __func__, __func__, eqp->name);

  if (n_previous < 0)
    n_previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? 1 : 0;

  _add_field(n_previous, eq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to all user-defined equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_user_create_fields(void)
{
  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];

    assert(eq != NULL);

    cs_equation_param_t  *eqp = eq->param;

    if (eqp->type !=  CS_EQUATION_TYPE_USER)
      continue;

    int  n_previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? 1 : 0;

    _add_field(n_previous, eq);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define the builder structure
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_define_builders(const cs_mesh_t       *mesh)
{
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Allocate and initialize a system builder */
    /* Not initialized here if it is a restart */

    if (eq->builder == NULL)
      eq->builder = cs_equation_builder_create(eq->param, mesh);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  }  /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define the context structure associated to each equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_define_context_structures(void)
{
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Allocate and initialize the context structures */
    /* Not initialized here if it is a restart */

    if (eq->scheme_context == NULL)
      eq->scheme_context = eq->init_context(eq->param,
                                            eq->field_id,
                                            eq->boundary_flux_id,
                                            eq->builder);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  }  /* Loop on equations */

  /* The system helper is defined during the definition of the scheme context.
   * When all system helper have been defined, one can allocate or update the
   * low-level assemble structures */

  cs_cdo_system_allocate_assembly();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a pointer to a core structure. If the input core structure is
 *         not allocated, then one allocates the structure.
 *
 * \param[in]       eq       pointer to a cs_equation_t structure
 * \param[in, out]  p_core   double pointer to a core structure to build
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_define_core(const cs_equation_t    *eq,
                        cs_equation_core_t    **p_core)
{
  cs_equation_core_t  *core = (p_core == NULL) ? NULL : *p_core;

  if (eq == NULL)
    return;

  if (core == NULL)
    BFT_MALLOC(core, 1, cs_equation_core_t);

  core->param = eq->param;
  core->builder = eq->builder;
  core->scheme_context = eq->scheme_context;

  *p_core = core;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initialize condition to all variable fields associated to
 *         each cs_equation_t structure.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_field_values(const cs_mesh_t             *mesh,
                              const cs_time_step_t        *ts)
{
  /* Loop on all equations */

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    const cs_equation_param_t  *eqp = eq->param;

#if defined(DEBUG) && !defined(NDEBUG)
    /* Check that the main structure for an equation have been defined */

    if (eq->builder == NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: A builder structure is expected for eq. \"%s\"\n",
                __func__, eqp->name);

    if (eq->scheme_context == NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: A context structure is expected for eq. \"%s\"\n",
                __func__, eqp->name);
#endif

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
 * \brief  Build the linear system for this equation (deprecated)
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

  /* Build the algebraic system to solve */

  eq->build_system(mesh, fld->val, eq->param, eq->builder, eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation (deprecated)
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
  cs_equation_builder_t  *eqb = eq->builder;
  cs_cdo_system_helper_t  *sh = eqb->system_helper;

  cs_real_t  *x = NULL;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  const cs_equation_param_t  *eqp = eq->param;
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_sles_t  *sles_param = eqp->sles_param;

  /* Handle parallelism (the the x array and for the rhs) */

  eq->prepare_solving(eq, &x);

  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);

  cs_sles_convergence_state_t code = cs_sles_solve(sles,
                                                   matrix,
                                                   sles_param->eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   sh->rhs,
                                                   x,
                                                   0,      /* aux. size */
                                                   NULL);  /* aux. buffers */

  if (sles_param->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "  <%s/sles_cvg> code %-d n_iters %d residual % -8.4e\n",
                  eqp->name, code, n_iters, residual);

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    const cs_range_set_t  *rset = cs_equation_get_range_set(eq);

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, // type and stride
                         x,
                         x);

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, // type and stride
                         sh->rhs,
                         sh->rhs);

  }

  /* Copy current field values to previous values */

  cs_field_current_to_previous(fld);

  /* Define the new field value for the current time */

  eq->update_field(x, sh->rhs, eq->param, eqb, eq->scheme_context, fld->val);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  /* Free memory */

  BFT_FREE(x);
  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);  /* free rhs and matrix */
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
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty equation structure", __func__);
  assert(cs_equation_uses_new_mechanism(eq));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
   * The linear solver is called inside and the field value is updated inside
   */

  eq->solve_steady_state(false, /* current to previous */
                         mesh,
                         eq->field_id,
                         eq->param,
                         eq->builder,
                         eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for an equation with an
 *         unsteady term
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in, out] eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(bool                        cur2prev,
                  const cs_mesh_t            *mesh,
                  cs_equation_t              *eq)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty equation structure", __func__);
  assert(cs_equation_uses_new_mechanism(eq));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
   * The linear solver is called inside and the field value is updated inside
   */

  eq->solve(cur2prev,
            mesh,
            eq->field_id,
            eq->param,
            eq->builder,
            eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for a steady-state equation.
 *         This is wrapper for the FORTRAN interface (limitation of the
 *         parameters to simple types).
 *
 * \param[in] eqname     name of the equation to solve
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_steady_state_wrapper(const char                 *eqname)
{
  cs_equation_t  *eq = cs_equation_by_name(eqname);

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty equation structure", __func__);
  assert(cs_equation_uses_new_mechanism(eq));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
     The linear solver is called inside and the field value is updated inside
  */

  eq->solve_steady_state(false, /* current to previous */
                         cs_glob_mesh,
                         eq->field_id,
                         eq->param,
                         eq->builder,
                         eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for an equation with an
 *         unsteady term. This is wrapper for the FORTRAN interface (limitation
 *         of the parameters to simple types)
 *
 * \param[in] cur2prev   true="current to previous" operation is performed
 * \param[in] eqname     name of the equation to solve
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_wrapper(bool                        cur2prev,
                          const char                 *eqname)
{
  cs_equation_t  *eq = cs_equation_by_name(eqname);

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty equation structure", __func__);
  assert(cs_equation_uses_new_mechanism(eq));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
     The linear solver is called inside and the field value is updated inside
  */
  eq->solve(cur2prev,
            cs_glob_mesh,
            eq->field_id,
            eq->param,
            eq->builder,
            eq->scheme_context);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the current to previous to all fields (and potentially arrays)
 *         related to an equation. This function fas to be called when a solve
 *         step is called with the parameter: cur2prev = false
 *
 * \param[in]   eq       pointer to a \ref cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_current_to_previous(const cs_equation_t    *eq)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty equation structure", __func__);

  if (eq->current_to_previous == NULL)
    return;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  eq->current_to_previous(eq->param, eq->builder, eq->scheme_context);

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
 *         cell of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of cell values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq,
                            bool                    previous)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *c_values = NULL;
  if (eq->get_cell_values != NULL)
    c_values = eq->get_cell_values(eq->scheme_context, previous);

  return c_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         face of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of face values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq,
                            bool                    previous)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *f_values = NULL;
  if (eq->get_face_values != NULL)
    f_values = eq->get_face_values(eq->scheme_context, previous);

  return f_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         edge of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of edge values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_edge_values(const cs_equation_t    *eq,
                            bool                    previous)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *e_values = NULL;
  if (eq->get_edge_values != NULL)
    e_values = eq->get_edge_values(eq->scheme_context, previous);

  return e_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         vertex of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of vertex values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_vertex_values(const cs_equation_t    *eq,
                              bool                    previous)
{
  if (eq == NULL)
    return NULL;

  cs_real_t  *v_values = NULL;
  if (eq->get_vertex_values != NULL)
    v_values = eq->get_vertex_values(eq->scheme_context, previous);

  return v_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over the domain of the current variable field
 *         associated to the given equation.
 *
 * \param[in]      connect    pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq         pointer to a \ref cs_equation_t structure
 * \param[in, out] integral   result of the computation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_integrate_variable(const cs_cdo_connect_t     *connect,
                               const cs_cdo_quantities_t  *cdoq,
                               const cs_equation_t        *eq,
                               cs_real_t                  *result)
{
  /* Initialize the value to return */

  *result = 0.;

  if (eq == NULL)
    return;

  const cs_equation_param_t  *eqp = eq->param;
  assert(eqp != NULL);
  if (eqp->dim > 1)
    bft_error(__FILE__, __LINE__, 0, "%s: (Eq. %s) Not implemented",
              __func__, eqp->name);

  /* Scalar-valued equation */

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *p_v = cs_equation_get_vertex_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

        cs_real_t  int_cell = 0.;
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          int_cell += cdoq->dcell_vol[j] * p_v[c2v->ids[j]];

        *result += int_cell;

      } /* Loop on cells */

    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *p_v = cs_equation_get_vertex_values(eq, false);
      const cs_real_t  *p_c = cs_equation_get_cell_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

        /* Shares between cell and vertex unknowns:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
         */

        cs_real_t  int_cell = 0.25*cdoq->cell_vol[c_id]*p_c[c_id];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          int_cell += 0.75*cdoq->dcell_vol[j] * p_v[c2v->ids[j]];

        *result += int_cell;

      } /* Loop on cells */

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:
    {
      const cs_real_t  *p_f = cs_equation_get_face_values(eq, false);
      const cs_real_t  *p_c = cs_equation_get_cell_values(eq, false);
      const cs_adjacency_t  *c2f = connect->c2f;

      assert(cdoq->pvol_fc != NULL); /* Sanity check */

      for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

        /* Shares between cell and face unknowns (assuming stab.coeff = 1/3):
           - the cell unknown stands for 1/4 of the cell volume
           - the face unknown stands for 3/4 of the dual cell volume
         */

        cs_real_t  int_cell = 0.25*cdoq->cell_vol[c_id]*p_c[c_id];
        for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
          int_cell += 0.75*cdoq->pvol_fc[j] * p_f[c2f->ids[j]];

        *result += int_cell;

      } /* Loop on cells */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: (Eq. %s). Not implemented.", __func__, eqp->name);

  } /* End of switch */

  /* Parallel synchronization */

  cs_parall_sum(1, CS_REAL_TYPE, result);
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
      const cs_real_t  *p_v = cs_equation_get_vertex_values(eq, false);

      cs_cdovb_scaleq_boundary_diff_flux(t_eval,
                                         eqp,
                                         p_v,
                                         eq->builder,
                                         eq->scheme_context,
                                         diff_flux);
    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *p_v = cs_equation_get_vertex_values(eq, false);
      const cs_real_t  *p_c = cs_equation_get_cell_values(eq, false);

      cs_cdovcb_scaleq_boundary_diff_flux(t_eval,
                                          eqp,
                                          p_v,
                                          p_c,
                                          eq->builder,
                                          eq->scheme_context,
                                          diff_flux);
    }
    break;

  case CS_SPACE_SCHEME_CDOFB:
    {
      const cs_real_t  *p_f = cs_equation_get_face_values(eq, false);
      const cs_real_t  *p_c = cs_equation_get_cell_values(eq, false);

      cs_cdofb_scaleq_boundary_diff_flux(t_eval,
                                         eqp,
                                         p_f,
                                         p_c,
                                         eq->builder,
                                         eq->scheme_context,
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

  const char  emsg[] = "%s: Computation of the diffusive and convective flux"
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
 * \brief Cellwise computation of the diffusive flux across the requested
 *        location. If the location is not the "natural" one (which depends on
 *        the space discretization scheme) then the diffusive flux is only an
 *        approximation.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      location    indicate where the flux has to be computed
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] diff_flux   value of the diffusive flux (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diffusive_flux(const cs_equation_t   *eq,
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

  case CS_SPACE_SCHEME_CDOFB:
    if (cs_flag_test(location, cs_flag_primal_face))
      cs_cdofb_scaleq_diff_flux_faces(fld->val,
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

  if (!(eqp->post_flag & CS_EQUATION_POST_PECLET))
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
  CS_UNUSED(connect);
  CS_UNUSED(cdoq);

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    assert(eq != NULL); /* Sanity check */
    const cs_equation_param_t  *eqp = eq->param;

    /* Check if the computation of the balance is requested */

    if (!(eqp->post_flag & CS_EQUATION_POST_BALANCE))
      continue;

    if (eq->compute_balance != NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Balance for equation %s is requested but\n"
                " this functionality is not available yet.\n",
                __func__, eqp->name);

    if (eq->main_ts_id > -1)    /* Activate timer statistics */
      cs_timer_stats_start(eq->main_ts_id);

    cs_cdo_balance_t  *b = eq->compute_balance(eqp,
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

    /* In case of postprocessing of the border faces, one has to check if there
       is a mesh modification. In particular, a removal of 2D extruded border
       faces*/

    bool  use_parent = (mesh->n_g_b_faces_all > mesh->n_g_b_faces) ?
      false : true;

    sprintf(postlabel, "%s.BdyFlux", eqp->name);

    /* Post-process the boundary fluxes (diffusive and convective) */

    cs_post_write_var(CS_POST_MESH_BOUNDARY,
                      CS_POST_WRITER_DEFAULT,
                      postlabel,
                      1,
                      true,                /* interlace */
                      use_parent,
                      CS_POST_TYPE_cs_real_t,
                      NULL,                /* values on cells */
                      NULL,                /* values at internal faces */
                      b->boundary_term,    /* values at border faces */
                      ts);                 /* time step structure */

    /* Free buffers */

    BFT_FREE(postlabel);
    cs_cdo_balance_destroy(&b);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cellwise stiffness matrix associated to the property
 *         given as a parameter and apply it to the pot array to define
 *         the resulting array associated to entities defined at loc_res
 *
 * \param[in]      eq        pointer to a \ref cs_equation_t structure
 * \param[in]      property  pointer to the property to consider
 * \param[in]      pot       array to multiply with the stiffness matrix
 * \param[in]      loc_res   location of entities in the resulting array
 * \param[in, out] res       resulting array
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_apply_stiffness(cs_equation_t          *eq,
                            const cs_property_t    *property,
                            const cs_real_t        *pot,
                            cs_flag_t               loc_res,
                            cs_real_t              *res)
{
  if (eq == NULL)
    return;

  /* Preliminary checkings */

  if (pot == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Input array not allocated.\n",
              __func__);
  if (res == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Resulting array not allocated.\n",
              __func__);
  if (eq->apply_stiffness == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Function not defined for this equation \"%s\".\n",
              __func__, cs_equation_get_name(eq));

  /* Perform the requested operation */

  if (eq->main_ts_id > -1)    /* Activate timer statistics */
    cs_timer_stats_start(eq->main_ts_id);

 eq->apply_stiffness(eq->param,
                      eq->builder,
                      eq->scheme_context,
                      property,
                      pot,
                      loc_res,
                      res);

 if (eq->main_ts_id > -1)
   cs_timer_stats_stop(eq->main_ts_id);
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

    eq->postprocess(eqp,
                    eq->builder,
                    eq->scheme_context);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
