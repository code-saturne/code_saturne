/*============================================================================
 * Routines to handle extended definitions of quantities (cs_xdef_t structures)
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
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_XDEF_DBG  0

/*============================================================================
 * Global static variables
 *============================================================================*/

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
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_volume_create(cs_xdef_type_t    type,
                      int               dim,
                      int               z_id,
                      cs_flag_t         state,
                      cs_flag_t         meta,
                      void             *input)
{
  cs_xdef_t  *d = NULL;

  BFT_MALLOC(d, 1, cs_xdef_t);

  d->type = type;
  d->support = CS_XDEF_SUPPORT_VOLUME;
  d->dim = dim;
  d->z_id = z_id;
  d->state = state;
  d->meta = meta;
  d->qtype = CS_QUADRATURE_BARY; // default value

  switch (type) {

  case CS_XDEF_BY_VALUE:
    {
      double  *_input = (double *)input;
      BFT_MALLOC(d->input, dim, double);
      double  *_input_cpy = (double *)d->input;

      for (int i = 0; i < dim; i++) _input_cpy[i] = _input[i];

      /* Update state flag */
      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_input_t  *a = (cs_xdef_analytic_input_t *)input;
      cs_xdef_analytic_input_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_analytic_input_t);
      b->func = a->func;
      b->input = a->input;

      d->input = b;
    }
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *a = (cs_xdef_array_input_t *)input;
      cs_xdef_array_input_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_array_input_t);
      b->stride = a->stride;
      b->loc = a->loc;
      b->values = a->values;
      b->is_owner = a->is_owner;
      b->index = a->index;

      /* Update state flag */
      if (cs_flag_test(b->loc, cs_flag_primal_cell) ||
          cs_flag_test(b->loc, cs_flag_dual_face_byc))
        d->state |= CS_FLAG_STATE_CELLWISE;

      d->input = b;
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)input;

      d->input = f;

      /* Update state flag */
      if (f->location_id == cs_mesh_location_get_id_by_name(N_("cells")))
        d->state |= CS_FLAG_STATE_CELLWISE;
    }
    break;

  case CS_XDEF_BY_QOV:
    {
      double  *_input = (double *)input;
      BFT_MALLOC(d->input, 1, double);
      double *_input_cpy = (double *)d->input;
      _input_cpy[0] = _input[0];
    }
    break;

  default:
    d->input = input;
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
 * \param[in]  input      pointer to a structure
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
                        void             *input)
{
  cs_xdef_t  *d = NULL;

  BFT_MALLOC(d, 1, cs_xdef_t);

  d->type = type;
  d->support = CS_XDEF_SUPPORT_BOUNDARY;
  d->dim = dim;
  d->z_id = z_id;
  d->state = state;
  d->meta = meta;
  d->qtype = CS_QUADRATURE_BARY; // default value

  switch (type) {

  case CS_XDEF_BY_VALUE:
    {
      double  *_input = (double *)input;
      BFT_MALLOC(d->input, dim, double);
      double  *_input_cpy = (double *)d->input;
      for (int i = 0; i < dim; i++) _input_cpy[i] = _input[i];

      /* Update state flag */
      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE;
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_input_t  *a = (cs_xdef_analytic_input_t *)input;
      cs_xdef_analytic_input_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_analytic_input_t);
      b->func = a->func;
      b->input = a->input;

      d->input = b;
    }
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *a = (cs_xdef_array_input_t *)input;
      cs_xdef_array_input_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_array_input_t);
      b->stride = a->stride;
      b->loc = a->loc;
      b->values = a->values;
      b->is_owner = a->is_owner;
      b->index = a->index;

      d->input = b;

      /* Update state flag */
      if (cs_flag_test(b->loc, cs_flag_primal_face))
        d->state |= CS_FLAG_STATE_FACEWISE;
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)input;

      d->input = &(f->id);
    }
    break;

  case CS_XDEF_BY_QOV:
    {
      double  *_input = (double *)input;
      BFT_MALLOC(d->input, 1, double);
      double  *_input_cpy = (double *)d->input;
      _input_cpy[0] = _input[0];

      /* Update state flag */
      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE;
    }
    break;

  default: // analytic functions or more generic functions
    d->input = input;
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
 * \param[in]  input      pointer to a structure storing the parameters (cast
 *                        on-the-fly according to the type of definition)
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_timestep_create(cs_xdef_type_t             type,
                        cs_flag_t                  state,
                        cs_flag_t                  meta,
                        void                      *input)
{
  cs_xdef_t  *d = NULL;

  BFT_MALLOC(d, 1, cs_xdef_t);

  d->type = type;
  d->support = CS_XDEF_SUPPORT_TIME;
  d->dim = 1;
  d->z_id = -1; // No associated zone
  d->state = state;
  d->meta = meta;
  d->qtype = CS_QUADRATURE_NONE; // default value

  switch (type) {

  case CS_XDEF_BY_VALUE:
    {
      double  *_input = (double *)input;
      BFT_MALLOC(d->input, 1, double);
      double  *_input_cpy = (double *)d->input;
      _input_cpy[0] = _input[0];

      /* Update state flag */
      d->state |= CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_STEADY;
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    {
      cs_xdef_timestep_input_t  *a = (cs_xdef_timestep_input_t *)input;
      cs_xdef_timestep_input_t  *b = NULL;

      BFT_MALLOC(b, 1, cs_xdef_timestep_input_t);
      b->func = a->func;
      b->input = a->input;

      d->input = b;
    }
    break;

  default:
    d->input = input;
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

  if (d->type == CS_XDEF_BY_ARRAY) {

    cs_xdef_array_input_t  *a = (cs_xdef_array_input_t *)d->input;
    if (a->is_owner)
      BFT_FREE(a->values);
    BFT_FREE(d->input);

  }
  else if (d->type == CS_XDEF_BY_TIME_FUNCTION ||
           d->type == CS_XDEF_BY_VALUE ||
           d->type == CS_XDEF_BY_ANALYTIC_FUNCTION ||
           d->type == CS_XDEF_BY_QOV)
    BFT_FREE(d->input);

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
                                src->input);
    break;

  case CS_XDEF_SUPPORT_TIME:
    cpy = cs_xdef_timestep_create(src->type,
                                  src->state,
                                  src->meta,
                                  src->input);
    break;

  case CS_XDEF_SUPPORT_BOUNDARY:
    cpy = cs_xdef_boundary_create(src->type,
                                  src->dim,
                                  src->z_id,
                                  src->state,
                                  src->meta,
                                  src->input);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case", __func__);

  }

  cpy->qtype = src->qtype;

  return cpy;
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

  cs_xdef_array_input_t  *a = (cs_xdef_array_input_t *)d->input;

  /* An array is already assigned and one manages the lifecycle */
  if (a->is_owner && a->values != NULL)
    BFT_FREE(a->values);

  /* Set the new values */
  a->is_owner = is_owner;
  a->values = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In case of definition by array, set the index to get access to the
 *         array values.
 *
 * \param[in, out]  d             pointer to a cs_xdef_t structure
 * \param[in]       array_index   index on array values
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_array_index(cs_xdef_t     *d,
                        cs_lnum_t     *array_index)
{
  if (d == NULL)
    return;

  if (d->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The given cs_xdef_t structure should be defined by array.",
              __func__);

  cs_xdef_array_input_t  *ai = (cs_xdef_array_input_t *)d->input;

  ai->index = array_index;
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
 * \brief  Output the settings related to a a cs_xdef_t structure
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log(const cs_xdef_t     *d)
{
  if (d == NULL)
    return;

  bool  is_uniform = false, is_steady = false, is_cellwise = false;
  if (d->state & CS_FLAG_STATE_UNIFORM)  is_uniform = true;
  if (d->state & CS_FLAG_STATE_STEADY)   is_steady = true;
  if (d->state & CS_FLAG_STATE_CELLWISE) is_cellwise = true;

  cs_log_printf(CS_LOG_SETUP,
                "  <Definition> uniform [%s], cellwise [%s], steady [%s],"
                " meta: %u\n",
                cs_base_strtf(is_uniform), cs_base_strtf(is_cellwise),
                cs_base_strtf(is_steady), d->meta);

  if (d->support == CS_XDEF_SUPPORT_VOLUME) {

    const cs_zone_t  *z = cs_volume_zone_by_id(d->z_id);
    assert(z != NULL);
    cs_log_printf(CS_LOG_SETUP,
                  "  <Definition> support: volume, zone: %d, %s,"
                  " mesh_location: %s\n",
                  z->id, z->name, cs_mesh_location_get_name(z->location_id));

  }
  else if (d->support == CS_XDEF_SUPPORT_BOUNDARY) {

    const cs_zone_t  *z = cs_boundary_zone_by_id(d->z_id);
    assert(z != NULL);
    cs_log_printf(CS_LOG_SETUP,
                  "  <Definition> support: boundary, zone: %d, %s,"
                  " mesh_location: %s\n",
                  z->id, z->name, cs_mesh_location_get_name(z->location_id));

  }
  else if (d->support == CS_XDEF_SUPPORT_TIME)
    cs_log_printf(CS_LOG_SETUP, " <Definition> support: time\n");

  switch (d->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_log_printf(CS_LOG_SETUP, "              by an analytical function\n");
    break;

  case CS_XDEF_BY_ARRAY:
    cs_log_printf(CS_LOG_SETUP, "              by an array\n");
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)d->input;

      cs_log_printf(CS_LOG_SETUP, "              by the field %s\n",
                    f->name);
    }
    break;

  case CS_XDEF_BY_FUNCTION:
    cs_log_printf(CS_LOG_SETUP, "              by function\n");
    break;

  case CS_XDEF_BY_QOV:
    cs_log_printf(CS_LOG_SETUP, "              by quantity over a volume\n");
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    cs_log_printf(CS_LOG_SETUP, "              by time function\n");
    break;

  case CS_XDEF_BY_VALUE:
    {
      cs_real_t *values = (cs_real_t *)d->input;

      if (d->dim == 1)
        cs_log_printf(CS_LOG_SETUP, "              by_value, % 5.3e\n",
                      values[0]);
      else if (d->dim == 3)
        cs_log_printf(CS_LOG_SETUP, "              by_value,"
                      " (% 5.3e, % 5.3e, % 5.3e)\n",
                      values[0], values[1], values[2]);
      else if (d->dim == 9)
        cs_log_printf(CS_LOG_SETUP, "              by_value,"
                      " ((% 4.2e, % 4.2e, % 4.2e) (% 4.2e, % 4.2e, % 4.2e)"
                      " (% 4.2e, % 4.2e, % 4.2e))\n",
                      values[0], values[1], values[2],
                      values[3], values[4], values[5],
                      values[6], values[7], values[8]);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid case. dim = %d (expected 3, 6 or 9)\n", d->dim);
    }
    break; // BY_VALUE


  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of description."));
    break;

  } /* switch on def_type */

  cs_log_printf(CS_LOG_SETUP, "  <Definition/Quadrature> %s\n",
                cs_quadrature_get_type_name(d->qtype));

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
