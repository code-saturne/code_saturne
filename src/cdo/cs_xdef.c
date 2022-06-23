/*============================================================================
 * Functions to handle extended definitions of quantities
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_field.h"
#include "cs_flag.h"
#include "cs_log.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

 /*!
   \file cs_xdef.c

   \brief Functions to handle extended definitions of quantities thanks to the
          cs_xdef_t structures.
 */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_XDEF_DBG  0

/*============================================================================
 * Global static variables
 *============================================================================*/

static const char *_xdef_type_name[]
  = {"CS_XDEF_BY_ANALYTIC_FUNCTION",
     "CS_XDEF_BY_ARRAY",
     "CS_XDEF_BY_DOF_FUNCTION",
     "CS_XDEF_BY_FIELD",
     "CS_XDEF_BY_FUNCTION",
     "CS_XDEF_BY_QOV",
     "CS_XDEF_BY_SUB_DEFINITIONS",
     "CS_XDEF_BY_TIME_FUNCTION",
     "CS_XDEF_BY_VALUE",
     "out of range"};

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on volumic
 *         elements
 *
 * \param[in]  type        type of definition
 * \param[in]  dim         dimension of the values to define
 * \param[in]  z_id        volume zone id
 * \param[in]  state       flag to know if this uniform, cellwise, steady...
 * \param[in]  meta        metadata associated to this description
 * \param[in]  context     pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_volume_create(cs_xdef_type_t           type,
                      int                      dim,
                      int                      z_id,
                      cs_flag_t                state,
                      cs_flag_t                meta,
                      void                    *context)
{
  cs_xdef_t  *d = NULL;

  BFT_MALLOC(d, 1, cs_xdef_t);

  d->type = type;
  d->support = CS_XDEF_SUPPORT_VOLUME;
  d->dim = dim;
  d->z_id = z_id;
  d->state = state;
  d->meta = meta;
  d->qtype = CS_QUADRATURE_BARY; /* default value */

  /* Now define the context pointer */

  switch (type) {

  case CS_XDEF_BY_VALUE:
    {
      double  *_context = (double *)context;
      BFT_MALLOC(d->context, dim, double);

      double  *_context_cpy = (double *)d->context;
      for (int i = 0; i < dim; i++) _context_cpy[i] = _context[i];

      /* Update state flag */
      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_context_t  *a = (cs_xdef_analytic_context_t *)context;
      cs_xdef_analytic_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_analytic_context_t);
      assert(a->z_id == z_id);
      b->z_id = a->z_id;
      b->func = a->func;
      b->input = a->input;
      b->free_input = a->free_input;

      d->context = b;
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    {
      cs_xdef_dof_context_t  *a = (cs_xdef_dof_context_t *)context;
      cs_xdef_dof_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_dof_context_t);
      assert(a->z_id == z_id);
      b->z_id = a->z_id;
      b->func = a->func;
      b->loc = a->loc;
      b->input = a->input;
      b->free_input = a->free_input;

      d->context = b;
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    {
      cs_xdef_time_func_context_t  *a = (cs_xdef_time_func_context_t *)context;
      cs_xdef_time_func_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_time_func_context_t);
      b->func = a->func;
      b->input = a->input;
      b->free_input = a->free_input;

      d->context = b;
    }
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_context_t  *a = (cs_xdef_array_context_t *)context;
      cs_xdef_array_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_array_context_t);
      assert(a->z_id == z_id);
      b->z_id = a->z_id;
      b->stride = a->stride;
      b->loc = a->loc;
      b->values = a->values;
      b->is_owner = a->is_owner;
      b->index = a->index;
      b->ids = a->ids;

      /* Update state flag */

      if (cs_flag_test(b->loc, cs_flag_primal_cell) ||
          cs_flag_test(b->loc, cs_flag_dual_face_byc))
        d->state |= CS_FLAG_STATE_CELLWISE;

      d->context = b;
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)context;

      d->context = f;
      assert(f != NULL);

      const cs_mesh_location_type_t  loc_type =
        cs_mesh_location_get_type(f->location_id);

      /* Update state flag */

      switch(loc_type) {

      case CS_MESH_LOCATION_CELLS:
        d->state |= CS_FLAG_STATE_CELLWISE;
        d->meta |= CS_FLAG_FULL_LOC;
        break;
      case CS_MESH_LOCATION_VERTICES:
        d->meta |= CS_FLAG_FULL_LOC;
        break;

      default:
        break; /* Nothing to do */
      }

    }
    break;

  case CS_XDEF_BY_QOV:
    {
      double  *_context = (double *)context;

      BFT_MALLOC(d->context, 1, double);

      double *_context_cpy = (double *)d->context;
      _context_cpy[0] = _context[0];
    }
    break;

  default: /* More generic functions e.g. CS_XDEF_BY_FUNCTION */
    d->context = context;  /* remark: context is used as an input structure.
                              The lifecycle of this pointer is not managed by
                              the current cs_xdef_t structure */
    break;
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on boundary
 *         elements
 *
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  context    pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_boundary_create(cs_xdef_type_t    type,
                        int               dim,
                        int               z_id,
                        cs_flag_t         state,
                        cs_flag_t         meta,
                        void             *context)
{
  cs_xdef_t  *d = NULL;

  BFT_MALLOC(d, 1, cs_xdef_t);

  d->type = type;
  d->support = CS_XDEF_SUPPORT_BOUNDARY;
  d->dim = dim;
  d->z_id = z_id;
  d->state = state;
  d->meta = meta;
  d->qtype = CS_QUADRATURE_BARY; /* default value */

  switch (type) {

  case CS_XDEF_BY_VALUE:
    {
      double  *_context = (double *)context;

      BFT_MALLOC(d->context, dim, double);

      double  *_context_cpy = (double *)d->context;
      for (int i = 0; i < dim; i++) _context_cpy[i] = _context[i];

      /* Update state flag */

      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE;
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_context_t  *a = (cs_xdef_analytic_context_t *)context;
      cs_xdef_analytic_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_analytic_context_t);
      assert(a->z_id == z_id);
      b->z_id = a->z_id;
      b->func = a->func;
      b->input = a->input;
      b->free_input = a->free_input;

      d->context = b;
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    {
      cs_xdef_dof_context_t  *a = (cs_xdef_dof_context_t *)context;
      cs_xdef_dof_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_dof_context_t);
      b->func = a->func;
      b->loc = a->loc;
      b->input = a->input;
      b->free_input = a->free_input;

      d->context = b;

      /* Update state flag */

      if (cs_flag_test(b->loc, cs_flag_primal_face))
        d->state |= CS_FLAG_STATE_FACEWISE;
    }
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_context_t  *a = (cs_xdef_array_context_t *)context;
      cs_xdef_array_context_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_array_context_t);
      b->z_id = a->z_id;
      b->stride = a->stride;
      b->loc = a->loc;
      b->values = a->values;
      b->is_owner = a->is_owner;
      b->index = a->index;
      b->ids = a->ids;

      d->context = b;

      /* Update state flag */

      if (cs_flag_test(b->loc, cs_flag_primal_face))
        d->state |= CS_FLAG_STATE_FACEWISE;
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)context;
      assert(f != NULL);
      d->context = f;

      const cs_mesh_location_type_t  loc_type =
        cs_mesh_location_get_type(f->location_id);

      /* Update flags */

      if (loc_type == CS_MESH_LOCATION_BOUNDARY_FACES) {
        d->meta |= CS_FLAG_FULL_LOC;
        d->state |= CS_FLAG_STATE_FACEWISE;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Definition by field on the boundary rely on a mesh"
                  " location defined at boundary faces.", __func__);
    }
    break;

  case CS_XDEF_BY_QOV:
    {
      double  *_context = (double *)context;

      BFT_MALLOC(d->context, 1, double);

      double  *_context_cpy = (double *)d->context;
      _context_cpy[0] = _context[0];

      /* Update state flag */

      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE;
    }
    break;

  default: /* More generic functions e.g. CS_XDEF_BY_FUNCTION */
    d->context = context;   /* remark: context is used as an input structure.
                               The lifecycle of this pointer is not managed by
                               the current cs_xdef_t structure */
    break;

  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure for setting the
 *         time step
 *
 * \param[in]  type       type of definition
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  context    pointer to a structure storing the parameters (cast
 *                        on-the-fly according to the type of definition)
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_timestep_create(cs_xdef_type_t       type,
                        cs_flag_t            state,
                        cs_flag_t            meta,
                        void                *context)
{
  cs_xdef_t  *d = NULL;

  BFT_MALLOC(d, 1, cs_xdef_t);

  d->type = type;
  d->support = CS_XDEF_SUPPORT_TIME;
  d->dim = 1;
  d->z_id = -1;                  /* no associated zone */
  d->state = state;
  d->meta = meta;
  d->qtype = CS_QUADRATURE_NONE; /* default value */

  switch (type) {

  case CS_XDEF_BY_VALUE:
    {
      double  *_context = (double *)context;

      BFT_MALLOC(d->context, 1, double);

      double  *_context_cpy = (double *)d->context;
      _context_cpy[0] = _context[0];

      /* Update state flag */

      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_STEADY;
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    {
      cs_xdef_time_func_context_t *a = (cs_xdef_time_func_context_t *)context;
      cs_xdef_time_func_context_t *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_time_func_context_t);
      b->func = a->func;
      b->input = a->input;
      b->free_input = a->free_input;

      d->state |= CS_FLAG_STATE_UNIFORM;
      d->context = b;
    }
    break;

  default:
    d->context = context;   /* remark: context is used as an input structure.
                               The lifecycle of this pointer is not managed by
                               the current cs_xdef_t structure */
    break;
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_xdef_t structure
 *
 * \param[in, out] d    pointer to a cs_xdef_t structure
 *
 * \return NULL
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_free(cs_xdef_t     *d)
{
  if (d == NULL)
    return d;

  switch (d->type) {

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_context_t  *a = (cs_xdef_array_context_t *)d->context;
      if (a->is_owner)
        BFT_FREE(a->values);

      /* ids and idx if set are only shared so that one does not have to free
         thme at this stage */

      BFT_FREE(d->context);
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_context_t *c = (cs_xdef_analytic_context_t *)d->context;

      if (c->free_input != NULL)
        c->input = c->free_input(c->input);

      BFT_FREE(d->context);
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    {
      cs_xdef_dof_context_t *c = (cs_xdef_dof_context_t *)d->context;

      if (c->free_input != NULL)
        c->input = c->free_input(c->input);

      BFT_FREE(d->context);
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    {
      cs_xdef_time_func_context_t *c =
        (cs_xdef_time_func_context_t *)d->context;

      if (c->free_input != NULL)
        c->input = c->free_input(c->input);

      BFT_FREE(d->context);
    }
    break;

  case CS_XDEF_BY_VALUE:
  case CS_XDEF_BY_QOV:
    BFT_FREE(d->context);
    break;

  default:
    break; /* Nothing special to do e.g. CS_XDEF_BY_FUNCTION */
  }

  BFT_FREE(d);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  copy a cs_xdef_t structure
 *
 * \param[in]  src    pointer to a cs_xdef_t structure to copy
 *
 * \return a pointer to a new allocated cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_copy(cs_xdef_t     *src)
{
  cs_xdef_t  *cpy = NULL;
  if (src == NULL)
    return cpy;

  /* In the case of a definition by array where the structure is not owner
     one sets the copy to be owner of the array in order to avoid a memory
     leak */

  switch (src->support) {

  case CS_XDEF_SUPPORT_VOLUME:
    cpy = cs_xdef_volume_create(src->type,
                                src->dim,
                                src->z_id,
                                src->state,
                                src->meta,
                                src->context);
    break;

  case CS_XDEF_SUPPORT_TIME:
    cpy = cs_xdef_timestep_create(src->type,
                                  src->state,
                                  src->meta,
                                  src->context);
    break;

  case CS_XDEF_SUPPORT_BOUNDARY:
    cpy = cs_xdef_boundary_create(src->type,
                                  src->dim,
                                  src->z_id,
                                  src->state,
                                  src->meta,
                                  src->context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case", __func__);

  }

  cpy->qtype = src->qtype;

  return cpy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief In the case of a definition by an analytic function, a time function
 *        or a function relying on degrees of freedom (DoFs), this function
 *        allows one to set a more or less complex input data structure.  This
 *        call should be done before the first evaluation call of the
 *        associated cs_xdef_t structure.
 *
 * \param[in, out]  d         pointer to a cs_xdef_t structure
 * \param[in]       input     pointer to an input structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_input_context(cs_xdef_t       *d,
                          void            *input)
{
  if (d == NULL)
    return;

  switch (d->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_context_t *c = (cs_xdef_analytic_context_t *)d->context;

      c->input = input;
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    {
      cs_xdef_dof_context_t *c = (cs_xdef_dof_context_t *)d->context;

      c->input = input;
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    {
      cs_xdef_time_func_context_t *c =
        (cs_xdef_time_func_context_t *)d->context;

      c->input = input;
    }
    break;

  default:
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  " %s: Setting a free input function is ignored.\n"
                  " The type of definition is not compatible.", __func__);
    break; /* Nothing special to do */

  } /* End of switch */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief In case of a definition by an analytic function, a time function or a
 *        function relying on degrees of freedom (DoFs). One can set a function
 *        to free a complex input data structure (please refer to \ref
 *        cs_xdef_free_input_t) for more details.
 *
 * \param[in, out]  d             pointer to a cs_xdef_t structure
 * \param[in]       free_input    pointer to a function which free the input
 *                                structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_free_input_function(cs_xdef_t               *d,
                                cs_xdef_free_input_t    *free_input)
{
  if (d == NULL)
    return;

  switch (d->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_context_t *c = (cs_xdef_analytic_context_t *)d->context;

      c->free_input = free_input;
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    {
      cs_xdef_dof_context_t *c = (cs_xdef_dof_context_t *)d->context;

      c->free_input = free_input;
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    {
      cs_xdef_time_func_context_t *c =
        (cs_xdef_time_func_context_t *)d->context;

      c->free_input = free_input;
    }
    break;

  default:
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  " %s: Setting a free input function is ignored.\n"
                  " The type of definition is not compatible.", __func__);
    break; /* Nothing special to do */

  } /* End of switch */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the array after having added
 *         this definition
 *
 * \param[in, out]  d          pointer to a cs_xdef_t structure
 * \param[in]       is_owner   manage or not the lifecycle of the array values
 * \param[in]       array      values
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array(cs_xdef_t     *d,
                  bool           is_owner,
                  cs_real_t     *array)
{
  if (d == NULL)
    return;

  if (d->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The given cs_xdef_t structure should be defined by array.",
              __func__);

  cs_xdef_array_context_t  *a = (cs_xdef_array_context_t *)d->context;

  /* An array is already assigned and one manages the lifecycle */

  if (a->is_owner && a->values != NULL)
    BFT_FREE(a->values);

  /* Set the new values */

  a->is_owner = is_owner;
  a->values = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the zone id related to the size
 *         of the array. By default, the zone id is the same as the zone id
 *         related to the definition so that there is no need to call this
 *         function.
 *
 * \param[in, out]  d       pointer to a cs_xdef_t structure
 * \param[in]       z_id    zone id associated to the array size
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array_zone_id(cs_xdef_t     *d,
                          int            z_id)
{
  if (d == NULL)
    return;

  if (d->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The given cs_xdef_t structure should be defined by array.",
              __func__);

  cs_xdef_array_context_t  *actx = d->context;

  actx->z_id = z_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the optional index and ids
 *         arrays that may be useful when operating on definitions by array
 *
 * \param[in, out]  d         pointer to a cs_xdef_t structure
 * \param[in]       index     optional pointer to an array of index values
 * \param[in]       ids       optional pointer to a list of entity ids
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array_pointers(cs_xdef_t            *d,
                           const cs_lnum_t      *index,
                           const cs_lnum_t      *ids)
{
  if (d == NULL)
    return;

  if (d->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The given cs_xdef_t structure should be defined by array.",
              __func__);

  cs_xdef_array_context_t  *actx = d->context;

  actx->index = index;
  actx->ids = ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for evaluating the given
 *         description
 *
 * \param[in, out]  d       pointer to a cs_xdef_t structure
 * \param[in]       qtype   type of quadrature
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_quadrature(cs_xdef_t              *d,
                       cs_quadrature_type_t    qtype)
{
  if (d == NULL)
    return;

  d->qtype = qtype;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the type of quadrature to use for evaluating the given
 *         description
 *
 * \param[in]  d       pointer to a cs_xdef_t structure
 *
 * \return the type of quadrature
 */
/*----------------------------------------------------------------------------*/

cs_quadrature_type_t
cs_xdef_get_quadrature(cs_xdef_t     *d)
{
  if (d == NULL)
    return CS_QUADRATURE_NONE;

  return  d->qtype;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the flag dedicated to the state
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 *
 * \return the value of the flag
 */
/*----------------------------------------------------------------------------*/

cs_xdef_type_t
cs_xdef_get_type(const cs_xdef_t     *d)
{
  if (d == NULL)
    return CS_N_XDEF_TYPES;
  else
    return d->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the flag dedicated to the state
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 *
 * \return the value of the flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_xdef_get_state_flag(const cs_xdef_t     *d)
{
  if (d == NULL)
    return 0;
  else
    return d->state;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a cs_xdef_t structure in the setup
 *         logging file
 *
 * \param[in] prefix    optional string
 * \param[in] d         pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log_setup(const char          *prefix,
                  const cs_xdef_t     *d)
{
  cs_xdef_log(CS_LOG_SETUP, prefix, d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a cs_xdef_t structure
 *
 * \param[in] log_type  related log file to consider
 * \param[in] prefix    optional string
 * \param[in] d         pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log(cs_log_t             log_type,
            const char          *prefix,
            const cs_xdef_t     *d)
{
  if (d == NULL)
    return;

  bool  is_uniform = false, is_steady = false, is_cellwise = false;
  if (d->state & CS_FLAG_STATE_UNIFORM)  is_uniform = true;
  if (d->state & CS_FLAG_STATE_STEADY)   is_steady = true;
  if (d->state & CS_FLAG_STATE_CELLWISE) is_cellwise = true;

  const char  *_p;
  const char _empty_prefix[2] = "";
  if (prefix == NULL)
    _p = _empty_prefix;
  else
    _p = prefix;

  cs_log_printf(log_type,
                "%s | Uniform %s Cellwise %s Steady %s Meta: %u\n",
                _p, cs_base_strtf(is_uniform), cs_base_strtf(is_cellwise),
                cs_base_strtf(is_steady), d->meta);

  /* Which support */
  /* ============= */

  if (d->support == CS_XDEF_SUPPORT_VOLUME) {

    const cs_zone_t  *z = cs_volume_zone_by_id(d->z_id);
    assert(z != NULL);
    cs_log_printf(log_type, "%s | Support:   volume | Zone: %s (id:%5d)\n",
                  _p, z->name, z->id);

  }
  else if (d->support == CS_XDEF_SUPPORT_BOUNDARY) {

    const cs_zone_t  *z = cs_boundary_zone_by_id(d->z_id);
    assert(z != NULL);
    cs_log_printf(log_type, "%s | Support: boundary | Zone: %s (id:%5d)\n",
                  _p, z->name, z->id);

  }
  else if (d->support == CS_XDEF_SUPPORT_TIME)
    cs_log_printf(log_type, "%s | Support: time\n", _p);

  /* Type of definition */
  /* ================== */

  switch (d->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_log_printf(log_type, "%s | Definition by an analytical function\n",
                  _p);
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    cs_log_printf(log_type, "%s | Definition by a DoF function\n", _p);
    break;

  case CS_XDEF_BY_ARRAY:
    cs_log_printf(log_type, "%s | Definition by an array\n", _p);
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)d->context;

      if (f == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  " Field pointer is set to NULL in a definition by field");

      cs_log_printf(log_type, "%s | Definition by the field \"%s\"\n",
                    _p, f->name);
    }
    break;

  case CS_XDEF_BY_FUNCTION:
    cs_log_printf(log_type, "%s | Definition by function\n", _p);
    break;

  case CS_XDEF_BY_QOV:
    cs_log_printf(log_type,
                  "%s | Definition by a quantity over a volume\n", _p);
    break;

  case CS_XDEF_BY_SUB_DEFINITIONS:
    cs_log_printf(log_type, "%s | Definition by sub-definitions\n", _p);
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    cs_log_printf(log_type, "%s | Definition by a time function\n", _p);
    break;

  case CS_XDEF_BY_VALUE:
    {
      cs_real_t *values = (cs_real_t *)d->context;

      if (d->dim == 1)
        cs_log_printf(log_type, "%s | Definition by_value: % 5.3e\n",
                      _p, values[0]);
      else if (d->dim == 3)
        cs_log_printf(log_type, "%s | Definition by_value:"
                      " [% 5.3e, % 5.3e, % 5.3e]\n",
                      _p, values[0], values[1], values[2]);
      else if (d->dim == 9)
        cs_log_printf(log_type, "%s | Definition by_value:"
                      " [[% 4.2e, % 4.2e, % 4.2e], [% 4.2e, % 4.2e, % 4.2e],"
                      " [% 4.2e, % 4.2e, % 4.2e]]\n",
                      _p, values[0], values[1], values[2], values[3], values[4],
                      values[5], values[6], values[7], values[8]);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid case. dim = %d (expected 3, 6 or 9)\n",
                  __func__, d->dim);
    }
    break; /* BY_VALUE */

  default:
    bft_error(__FILE__, __LINE__, 0, _("%s: Invalid type of description."),
              __func__);
    break;

  } /* switch on def_type */

  cs_log_printf(log_type, "%s | Quadrature: %s\n",
                _p, cs_quadrature_get_type_name(d->qtype));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the cs_xdef_type's name string
 *
 * \param[in] xdef_type  type to query
 *
 * \return a pointer to mathing name string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_xdef_type_get_name(cs_xdef_type_t  xdef_type)
{
  if (xdef_type < 0 || xdef_type >= CS_N_XDEF_TYPES)
    xdef_type = CS_N_XDEF_TYPES;

  return _xdef_type_name[xdef_type];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
