/*============================================================================
 * Functions and structures to deal with source term computations
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh_location.h"
#include "cs_evaluate.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_SOURCE_TERM_DBG 0

struct _cs_source_term_t {

  char *restrict name;      // short description of the source term
  int            ml_id;     // id of the related mesh location structure

  /* Specification related to the way of computing the source term */
  cs_source_term_reduction_t  reduc;      // where reduction is done
  cs_source_term_type_t       type;       // mass, head loss...
  cs_quadra_type_t            quad_type;  // barycentric, higher, highest
  cs_param_var_type_t         var_type;   // scalar, vector...
  cs_param_def_type_t         def_type;   // by value, by function...
  cs_def_t                    def;

  /* Useful buffers to deal with more complex definitions */
  cs_desc_t    array_desc;  // short description of the related array
  cs_real_t   *array;       // if the source term hinges on an array
  const void  *struc;       // if one needs a structure. Only shared

};

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_st[] =
  " Stop setting an empty cs_source_term_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                   const cs_cdo_connect_t       *connect,
                                   const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
  cs_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a cs_source_term_t structure
 *
 * \param[in] st_name     name of the related source term
 * \param[in] ml_id       id of the related mesh location
 * \param[in] st_type     type of source term to create
 * \param[in] red_type    type of reduction to apply
 * \param[in] var_type    type of variables (scalar, vector, tensor...)
 *
 * \return a pointer to a new allocated source term structure
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_source_term_create(const char                  *name,
                      int                          ml_id,
                      cs_source_term_type_t        st_type,
                      cs_source_term_reduction_t   red_type,
                      cs_param_var_type_t          var_type)
{
  cs_source_term_t  *st = NULL;

  BFT_MALLOC(st, 1, cs_source_term_t);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(st->name, len, char);
  strncpy(st->name, name, len);

  st->type = st_type;
  st->var_type = var_type;
  st->ml_id = ml_id;

  /* Default initialization */
  st->reduc = red_type;
  st->def_type = CS_PARAM_N_DEF_TYPES;
  st->quad_type = CS_QUADRATURE_BARY;
  st->def.get.val = 0.;

  st->array_desc.location = 0;
  st->array_desc.state = 0;
  st->array = NULL;

  // st->struc is only shared when used

  return st;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_source_term_free(cs_source_term_t   *st)
{
  if (st == NULL)
    return st;

  BFT_FREE(st->name);

  if (st->array_desc.state & CS_FLAG_STATE_OWNER) {
    if (st->array != NULL)
      BFT_FREE(st->array);
  }

  BFT_FREE(st);

  return NULL;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set advanced parameters which are defined by default in a
 *         source term structure.
 *
 * \param[in, out]  st          pointer to a cs_source_term_t structure
 * \param[in]       quad_type   type of quadrature to use
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_quadrature(cs_source_term_t  *st,
                              cs_quadra_type_t   quad_type)
{
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  st->quad_type = quad_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set advanced parameters which are defined by default in a
 *         source term structure.
 *
 * \param[in, out]  st        pointer to a cs_source_term_t structure
 * \param[in]       red_type  type of reduction to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_reduction(cs_source_term_t             *st,
                             cs_source_term_reduction_t    red_type)
{
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  st->reduc = red_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the name related to a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return the name of the source term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_source_term_get_name(const cs_source_term_t   *st)
{
  if (st == NULL)
    return NULL;

  return st->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the type of reduction applied to a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return the type of reduction
 */
/*----------------------------------------------------------------------------*/

cs_source_term_reduction_t
cs_source_term_get_reduction(const cs_source_term_t   *st)
{
  if (st == NULL)
    return CS_N_SOURCE_TERM_REDUCTIONS;

  return st->reduc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the content of a cs_source_term_t structure
 *
 * \param[in] eqname  name of the related equation
 * \param[in] st      pointer to a cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_summary(const char               *eqname,
                       const cs_source_term_t   *st)
{
  const char  *eqn, _eqn[] = "Equation";

  if (eqname == NULL)
    eqn = _eqn;
  else
    eqn = eqname;

  if (st == NULL) {
    bft_printf("  <%s/NULL>\n", eqn);
    return;
  }

  bft_printf("  <%s/%s> type: ", eqn, st->name);
  switch (st->type) {

  case CS_SOURCE_TERM_GRAVITY:
    bft_printf(" Gravity");
    break;

  case CS_SOURCE_TERM_USER:
    bft_printf(" User-defined");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid type of source term.");
    break;

  }
  bft_printf(" mesh_location: %s\n", cs_mesh_location_get_name(st->ml_id));

  bft_printf("  <%s/%s> Definition: %s\n",
             eqn, st->name, cs_param_get_def_type_name(st->def_type));
  if (st->def_type == CS_PARAM_DEF_BY_ANALYTIC_FUNCTION)
    bft_printf("  <%s/%s> Quadrature: %s\n",
               eqn, st->name, cs_quadrature_get_type_name(st->quad_type));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic way to define the value of a cs_source_term_t structure
 *
 * \param[in, out]  pty     pointer to a cs_source_term_t structure
 * \param[in]       val     accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_value(cs_source_term_t    *st,
                            const char          *val)
{
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  st->def_type = CS_PARAM_DEF_BY_VALUE;

  switch (st->var_type) {

  case CS_PARAM_VAR_SCAL:
    cs_param_set_get(CS_PARAM_VAR_SCAL, (const void *)val, &(st->def.get));
    break;

  case CS_PARAM_VAR_VECT:
    cs_param_set_get(CS_PARAM_VAR_VECT, (const void *)val, &(st->def.get));
    break;

  case CS_PARAM_VAR_TENS:
    cs_param_set_get(CS_PARAM_VAR_TENS, (const void *)val, &(st->def.get));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Invalid type of source term."));
    break;

  } /* switch on variable type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_source_term_t structure thanks to an analytic function
 *
 * \param[in, out]  st      pointer to a cs_source_term_t structure
 * \param[in]       func    pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_analytic(cs_source_term_t      *st,
                               cs_analytic_func_t    *func)
{
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  st->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  st->def.analytic = func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_source_term_t structure thanks to an array of values
 *
 * \param[in, out]  st       pointer to a cs_source_term_t structure
 * \param[in]       desc     description of the main feature of this array
 * \param[in]       array    pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_array(cs_source_term_t    *st,
                            cs_desc_t            desc,
                            cs_real_t           *array)
{
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  st->def_type = CS_PARAM_DEF_BY_ARRAY;
  st->array_desc.location = desc.location;
  st->array_desc.state = desc.state;
  st->array = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term
 *
 * \param[in]      dof_desc   description of the associated DoF
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute(cs_desc_t                     dof_desc,
                       const cs_source_term_t       *source,
                       double                       *p_values[])
{
  const int  stride = 1; // Only this case is managed up to now
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  double  *values = *p_values;

  if (source == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  cs_lnum_t n_ent = 0;
  if (cs_cdo_same_support(dof_desc.location, cs_cdo_dual_cell) ||
      cs_cdo_same_support(dof_desc.location, cs_cdo_primal_vtx))
    n_ent = quant->n_vertices;
  else if (cs_cdo_same_support(dof_desc.location, cs_cdo_primal_cell))
    n_ent = quant->n_cells;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid case. Not able to compute an evaluation.\n"));

  /* Initialize values */
  if (values == NULL)
    BFT_MALLOC(values, n_ent*stride, double);
  for (cs_lnum_t i = 0; i < n_ent*stride; i++)
    values[i] = 0.0;

  if (dof_desc.state & CS_FLAG_STATE_POTENTIAL) {

    switch (source->def_type) {

    case CS_PARAM_DEF_BY_VALUE:
      cs_evaluate_potential_from_value(dof_desc.location,
                                       source->ml_id,
                                       source->def.get,
                                       values);
      break;

    case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
      cs_evaluate_potential_from_analytic(dof_desc.location,
                                          source->ml_id,
                                          source->def.analytic,
                                          values);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, _(" Invalid type of definition.\n"));

    } /* Switch according to def_type */

  }
  else if (dof_desc.state & CS_FLAG_STATE_DENSITY) {

    switch (source->def_type) {

    case CS_PARAM_DEF_BY_VALUE:
      cs_evaluate_density_from_value(dof_desc.location,
                                     source->ml_id,
                                     source->def.get,
                                     values);
      break;

    case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
      cs_evaluate_density_from_analytic(dof_desc.location,
                                        source->ml_id,
                                        source->def.analytic,
                                        source->quad_type,
                                        values);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, _(" Invalid type of definition.\n"));

    } /* Switch according to def_type */

  } /* Density variable */

  /* Return values */
  *p_values = values;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
