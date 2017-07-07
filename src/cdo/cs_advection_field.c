/*============================================================================
 * Manage the definition/setting of advection fields
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_xdef.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

#define CS_ADVECTION_FIELD_DBG  0

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_adv[] =
  " Stop setting an empty cs_adv_field_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

  /* Advection fields attached to the computational domain */
static int  _n_adv_fields = 0;
static cs_adv_field_t  **_adv_fields = NULL;


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
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_shared_pointers(const cs_cdo_quantities_t  *quant,
                                       const cs_cdo_connect_t     *connect,
                                       const cs_time_step_t       *time_step)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
  cs_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of allocated cs_adv_field_t structures
 *
 * \return the number of advection fields
 */
/*----------------------------------------------------------------------------*/

int
cs_advection_field_get_n_fields(void)
{
  return _n_adv_fields;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search in the array of advection field structures which one has
 *         the name given in argument
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to a cs_adv_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_by_name(const char   *name)
{
  if (_n_adv_fields <= 0)
    return NULL;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];
    if (cs_advection_field_check_name(adv, name))
      return adv;

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search in the array of advection field structures which one has
 *         the id given in argument
 *
 * \param[in]  id        identification number
 *
 * \return a pointer to a cs_adv_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_by_id(int      id)
{
  if (_n_adv_fields <= 0)
    return NULL;
  if (id >= _n_adv_fields || id < 0)
    return NULL;
  if (_adv_fields == NULL)
    return NULL;

  return _adv_fields[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add and initialize a new advection field structure
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to the new allocated cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_add(const char   *name)
{
  if (name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " A non-empty name is mandatory to add a new advection field");

  cs_adv_field_t  *adv = cs_advection_field_by_name(name);
  if (adv != NULL) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" An existing advection field has already the name %s.\n"
                    " Stop adding this advection field.\n"), name);
    return adv;
  }

  int  new_id = _n_adv_fields;
  _n_adv_fields++;
  BFT_REALLOC(_adv_fields, _n_adv_fields, cs_adv_field_t *);
  _adv_fields[new_id] = NULL;

  BFT_MALLOC(adv, 1, cs_adv_field_t);

  adv->id = new_id;

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(adv->name, len, char);
  strncpy(adv->name, name, len);

  /* Default initialization */
  adv->loc_flag = 0;
  adv->flag = 0;
  adv->vtx_field_id = -1;
  adv->cell_field_id = -1;
  adv->definition = NULL;

  /* Function pointers */
  adv->get_eval_all_vertices = NULL;
  adv->get_eval_at_cell = NULL;
  adv->get_eval_at_cell_cw = NULL;
  adv->get_eval_at_xyz_cw = NULL;

  /* Store the new advection field */
  _adv_fields[new_id] = adv;

  return adv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all alllocated cs_adv_field_t structures and its related array
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_destroy_all(void)
{
  if (_adv_fields == NULL)
    return;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];

    adv->definition = cs_xdef_free(adv->definition);

    BFT_FREE(adv->name);
    BFT_FREE(adv);

    /* All other pointers are shared */

  } // Loop on advection fields

  BFT_FREE(_adv_fields);
  _n_adv_fields = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given advection field has the name ref_name
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure to test
 * \param[in]  ref_name    name of the advection field to find
 *
 * \return true if the name of the advection field is ref_name otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_check_name(const cs_adv_field_t   *adv,
                              const char             *ref_name)
{
  if (adv == NULL)
    return false;

  int  reflen = strlen(ref_name);
  int  len = strlen(adv->name);

  if (reflen == len) {
    if (strcmp(ref_name, adv->name) == 0)
      return true;
    else
      return false;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the advection field is uniform, otherwise false
 *
 * \param[in]    adv    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_is_uniform(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return false;

  if (adv->definition->state & CS_FLAG_STATE_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the advection field is uniform in each cell
 *         otherwise false
 *
 * \param[in]    adv    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_is_cellwise(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return false;

  cs_flag_t  state = adv->definition->state;

  if (state & CS_FLAG_STATE_UNIFORM)
    return true;
  if (state & CS_FLAG_STATE_CELLWISE)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of an advection field
 *
 * \param[in]    adv    pointer to an advection field structure
 *
 * \return  the name of the related advection field
 */
/*----------------------------------------------------------------------------*/

const char *
cs_advection_field_get_name(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return NULL;

  return adv->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of definition used to set the current advection
 *         field structure
 *
 * \param[in]    adv    pointer to an advection field structure
 *
 * \return  the type of definition
 */
/*----------------------------------------------------------------------------*/

cs_xdef_type_t
cs_advection_field_get_deftype(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return CS_N_XDEF_TYPES;

  return cs_xdef_get_type(adv->definition);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print all setup information related to cs_adv_field_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_log_setup(void)
{
  if (_adv_fields == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the advection field\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_advection_fields       %d\n",
                _n_adv_fields);

  for (int i = 0; i < _n_adv_fields; i++) {

    const cs_adv_field_t  *adv = _adv_fields[i];

    cs_log_printf(CS_LOG_SETUP, " <AdvectionField/%s> id: %d\n",
                  adv->name, adv->id);
    cs_xdef_log(adv->definition);

  } // Loop on advection fields

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       key       key related to the member of adv to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_option(cs_adv_field_t            *adv,
                              cs_advection_field_key_t   key)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  switch(key) {

  case CS_ADVKEY_DEFINE_AT_CELLS:
    adv->loc_flag |= CS_FLAG_CELL;
    break;
  case CS_ADVKEY_DEFINE_AT_VERTICES:
    adv->loc_flag |= CS_FLAG_VERTEX;
    break;
  case CS_ADVKEY_POST_COURANT:
    adv->flag |= CS_ADVECTION_FIELD_POST_COURANT;
    break;
  case CS_ADVKEY_STATE_STEADY:
    adv->flag |= CS_ADVECTION_FIELD_STEADY;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key not implemented for setting an advection field."));
    break;

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       vector    values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_value(cs_adv_field_t    *adv,
                                cs_real_t          vector[3])
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0;

  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                          3, // dim.
                                          0, // zone_id = 0 => all cells
                                          state_flag,
                                          meta_flag,
                                          vector);

  /* Set function pointers */
  adv->get_eval_all_vertices = cs_xdef_eval_vector_by_val;
  adv->get_eval_at_cell = cs_xdef_eval_vector_by_val;
  adv->get_eval_at_cell_cw = cs_xdef_eval_cw_vector_by_val;
  adv->get_eval_at_xyz_cw = cs_xdef_eval_cw_vector_at_xyz_by_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an analytic function
 *
 * \param[in, out]  adv     pointer to a cs_adv_field_t structure
 * \param[in]       func    pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_analytic(cs_adv_field_t        *adv,
                                   cs_analytic_func_t    *func)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0;

  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          3, // dim.
                                          0, // zone_id = 0 => all cells
                                          state_flag,
                                          meta_flag,
                                          (void *)func);

  /* Set function pointers */
  adv->get_eval_all_vertices = cs_xdef_eval_at_vertices_by_analytic;
  adv->get_eval_at_cell = cs_xdef_eval_at_cells_by_analytic;
  adv->get_eval_at_cell_cw = cs_xdef_eval_cw_cell_by_analytic;
  adv->get_eval_at_xyz_cw = cs_xdef_eval_cw_at_xyz_by_analytic;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_array(cs_adv_field_t    *adv,
                                cs_flag_t          loc,
                                cs_real_t         *array,
                                cs_lnum_t         *index)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0; // Will be updated during the creation
  cs_flag_t  meta_flag = 0;  // metadata
  cs_xdef_array_input_t  input = {.stride = 3,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };

  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                          3, // dim
                                          0, // zone_id
                                          state_flag,
                                          meta_flag,
                                          &input);

  /* Set function pointers */
  adv->get_eval_all_vertices = cs_xdef_eval_3_at_all_vertices_by_array;
  adv->get_eval_at_cell = cs_xdef_eval_nd_at_cells_by_array;
  adv->get_eval_at_cell_cw = cs_xdef_eval_cw_cell_by_array;
  adv->get_eval_at_xyz_cw = cs_xdef_eval_cw_3_at_xyz_by_array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_field(cs_adv_field_t    *adv,
                                cs_field_t        *field)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0; // Will be updated during the creation
  cs_flag_t  meta_flag = 0;  // metadata

  assert(field->dim == 3); /* sanity check since fields are either at vertices
                              or at cells in this case */
  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_FIELD,
                                          3,
                                          0, // zone_id
                                          state_flag,
                                          meta_flag,
                                          field);

  /* Set function pointers */
  adv->get_eval_all_vertices = NULL;
  adv->get_eval_at_cell = cs_xdef_eval_cell_by_field;
  adv->get_eval_at_cell_cw = cs_xdef_eval_cw_cell_by_field;
  adv->get_eval_at_xyz_cw = cs_xdef_eval_cw_3_at_xyz_by_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create all needed cs_field_t structures related to an advection
 *         field
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_create_fields(void)
{
  int  len;
  char *field_name = NULL;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];
    assert(adv != NULL);

    bool  has_previous = (adv->flag & CS_ADVECTION_FIELD_STEADY) ? true:false;
    int  field_mask = CS_FIELD_PROPERTY;

    if (adv->loc_flag & CS_FLAG_VERTEX) { // Add a field attached to vertices

      /* Define the name of the field */
      len = strlen(adv->name) + strlen("_vertices") + 1;
      BFT_MALLOC(field_name, len, char);
      sprintf(field_name, "%s_vertices", adv->name);

      cs_field_t  *fld = cs_field_create(field_name,
                                         field_mask,
                                         CS_MESH_LOCATION_VERTICES,
                                         3,    // always a vector-valued field
                                         has_previous);

      cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

      adv->vtx_field_id = cs_field_id_by_name(field_name);

      BFT_FREE(field_name);

    } // Add a field attached to vertices

    if (adv->loc_flag & CS_FLAG_CELL) { // Add a field attached to cells

      /* Define the name of the field */
      len = strlen(adv->name) + strlen("_cells") + 1;
      BFT_MALLOC(field_name, len, char);
      sprintf(field_name, "%s_cells", adv->name);

      cs_field_t  *fld = cs_field_create(field_name,
                                         field_mask,
                                         CS_MESH_LOCATION_CELLS,
                                         3,    // always a vector-valued field
                                         has_previous);

      cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

      adv->cell_field_id = cs_field_id_by_name(field_name);

      BFT_FREE(field_name);

    } // Add a field attached to vertices

  } // Loop on advection fields
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a cs_field_t structure related to an advection field and a mesh
 *         location
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure
 * \param[in]  ml_type     type of mesh location (cells or vertices)
 *
 * \return a pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_advection_field_get_field(cs_adv_field_t           *adv,
                             cs_mesh_location_type_t   ml_type)
{
  cs_field_t  *f = NULL;
  if (adv == NULL)
    return f;

  switch (ml_type) {
  case CS_MESH_LOCATION_CELLS:
    assert(adv->cell_field_id > -1);
    f = cs_field_by_id(adv->cell_field_id);
    break;

  case CS_MESH_LOCATION_VERTICES:
    assert(adv->cell_field_id > -1);
    f = cs_field_by_id(adv->vtx_field_id);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid mesh location type to retrieve an advection field.\n");
  }

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at the cell center
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_adv_field_t   *adv,
                           cs_nvec3_t             *vect)
{
  /* Initialize the vector */
  vect->meas = 0.;
  for (int k = 0; k < 3; k++)
    vect->unitv[k] = 0;

  if (adv == NULL)
    return;

  cs_real_3_t  vector_values = {0, 0, 0};
  cs_xdef_t  *def = adv->definition;

  adv->get_eval_at_cell_cw(cm, cs_time_step, def->input, vector_values);

  cs_nvec3(vector_values, vect);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at the cell center
 *
 * \param[in]      c_id    id of the current cell
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_cell_vector(cs_lnum_t               c_id,
                                   const cs_adv_field_t   *adv,
                                   cs_nvec3_t             *vect)
{
  /* Initialize the vector */
  vect->meas = 0.;
  for (int k = 0; k < 3; k++)
    vect->unitv[k] = 0;

  if (adv == NULL)
    return;

  cs_real_3_t  vector_values = {0, 0, 0};
  cs_xdef_t  *def = adv->definition;

  adv->get_eval_at_cell(1,
                        &c_id,
                        true,
                        cs_glob_mesh,
                        cs_cdo_connect,
                        cs_cdo_quant,
                        cs_time_step,
                        def->input,
                        vector_values);

  cs_nvec3(vector_values, vect);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field for a given face
 *
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      xyz     coordinates where to evaluate the advection field
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_at_xyz(const cs_adv_field_t   *adv,
                              const cs_cell_mesh_t   *cm,
                              const cs_real_3_t       xyz,
                              cs_nvec3_t             *vect)
{
  /* Initialize the vector */
  vect->meas = 0.;
  for (int k = 0; k < 3; k++)
    vect->unitv[k] = 0;

  if (adv == NULL)
    return;

  cs_real_3_t  vector_values = {0, 0, 0};
  cs_xdef_t  *def = adv->definition;

  /* Sanity checks */
  assert(def != NULL);
  if (adv->get_eval_at_xyz_cw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Can not evaluate an advection field at xyz.");

  adv->get_eval_at_xyz_cw(cm, 1, xyz, cs_time_step, def->input,
                          vector_values);

  cs_nvec3(vector_values, vect);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at cell centers
 *
 * \param[in]      adv           pointer to a cs_adv_field_t structure
 * \param[in, out] cell_values   array of values at cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_cells(const cs_adv_field_t  *adv,
                            cs_real_t             *cell_values)
{
  if (adv == NULL)
    return;

  cs_xdef_t  *def = adv->definition;

  /* Sanity checks */
  assert(def != NULL);
  if (adv->get_eval_at_cell == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Can not evaluate the advection field.", __func__);

  adv->get_eval_at_cell(cs_cdo_quant->n_cells,
                        NULL,
                        false, // no compact
                        cs_glob_mesh,
                        cs_cdo_connect,
                        cs_cdo_quant,
                        cs_time_step,
                        def->input,
                        cell_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at vertices
 *
 * \param[in]      adv          pointer to a cs_adv_field_t structure
 * \param[in, out] vtx_values   array storing the results
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_vertices(const cs_adv_field_t  *adv,
                               cs_real_t             *vtx_values)
{
  if (adv == NULL)
    return;

  cs_xdef_t  *def = adv->definition;

  /* Sanity checks */
  assert(def != NULL);
  if (adv->get_eval_all_vertices == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Can not evaluate the advection field.", __func__);

  adv->get_eval_all_vertices(cs_cdo_quant->n_vertices,
                             NULL,
                             false, // no compact
                             cs_glob_mesh,
                             cs_cdo_connect,
                             cs_cdo_quant,
                             cs_time_step,
                             def->input,
                             vtx_values);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         the dual faces of a cell
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      adv      pointer to a cs_adv_field_t structure
 * \param[in, out] fluxes   array of values attached to dual faces of a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_flux_dfaces(const cs_cell_mesh_t         *cm,
                                   const cs_adv_field_t         *adv,
                                   cs_real_t                    *fluxes)
{
  if (adv == NULL)
    return;

  if (fluxes == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " fluxes array should be allocated before the call.");

  cs_xdef_t  *def = adv->definition;

  /* Sanity checks */
  assert(def != NULL);

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      /* Retrieve the advection field: Switch to a cs_nvec3_t representation */
      cs_real_3_t  cell_vector;
      cs_xdef_eval_cw_vector_by_val(cm, cs_time_step, def->input, cell_vector);
      cs_nvec3_t  adv_vect;
      cs_nvec3(cell_vector, &adv_vect);

      /* Sanity check */
      assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PEQ));

      /* Loop on cell edges */
      for (short int e = 0; e < cm->n_ec; e++)
        fluxes[e] = adv_vect.meas * cm->dface[e].meas *
          _dp3(adv_vect.unitv, cm->dface[e].unitv);

    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_EFQ |
                          CS_CDO_LOCAL_PFQ));

      cs_quadrature_type_t  qtype = cs_xdef_get_quadrature(def);

      /* Loop on cell edges */
      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_quant_t  edge = cm->edge[e];

        // Two triangles composing the dual face inside a cell
        const short int  f0 = cm->e2f_ids[2*e];
        const cs_nvec3_t  sef0 = cm->sefc[2*e];
        const cs_quant_t  qf0 = cm->face[f0];

        const short int  f1 = cm->e2f_ids[2*e+1];
        const cs_nvec3_t  sef1 = cm->sefc[2*e+1];
        const cs_quant_t  qf1 = cm->face[f1];

        fluxes[e] = 0.;
        switch (qtype) {

        case CS_QUADRATURE_NONE:
        case CS_QUADRATURE_BARY:
        case CS_QUADRATURE_BARY_SUBDIV:
          {
            cs_real_3_t  xg[2], adv_xg[2];

            for (int k = 0; k < 3; k++) {
              const double  xec = cm->xc[k] + edge.center[k];
              xg[0][k] = xec + qf0.center[k];
              xg[0][k] *= cs_math_onethird;
              xg[1][k] = xec + qf1.center[k];
              xg[1][k] *= cs_math_onethird;
            }

            cs_xdef_eval_cw_at_xyz_by_analytic(cm,
                                               2, (const cs_real_t *)xg,
                                               cs_time_step,
                                               def->input,
                                               (cs_real_t *)adv_xg);

            fluxes[e] = sef0.meas * _dp3(adv_xg[0], sef0.unitv)
              + sef1.meas * _dp3(adv_xg[1], sef1.unitv);
          }
          break;

        case CS_QUADRATURE_HIGHER:
          {
            cs_real_t  w[2];
            cs_real_3_t  gpts[6], eval[6];

            // Two triangles composing the dual face inside a cell
            cs_quadrature_tria_3pts(edge.center, qf0.center, cm->xc,
                                    sef0.meas,
                                    gpts, w);

            /* Evaluate the field at the three quadrature points */
            cs_quadrature_tria_3pts(edge.center, qf1.center, cm->xc,
                                    sef1.meas,
                                    gpts + 3, w + 1);

            cs_xdef_eval_cw_at_xyz_by_analytic(cm,
                                               6, (const cs_real_t *)gpts,
                                               cs_time_step,
                                               def->input,
                                               (cs_real_t *)eval);

            cs_real_t  add0 = 0, add1 = 0;
            for (int p = 0; p < 3; p++) add0 += _dp3(eval[p], sef0.unitv);
            add0 *= w[0];
            for (int p = 0; p < 3; p++) add1 += _dp3(eval[p+3], sef1.unitv);
            add1 *= w[1];

            fluxes[e] = add0 + add1;
          }
          break;

        case CS_QUADRATURE_HIGHEST: // Not yet implemented
        default:
          bft_error(__FILE__, __LINE__, 0, " Invalid type of quadrature.");
          break;

        } // switch type of quadrature

      } // Loop on cell edges

    } // definition by analytic function
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *input = (cs_xdef_array_input_t *)def->input;

      /* Test if location has at least the pattern of the reference support */
      if (cs_test_flag(input->loc, cs_cdo_dual_face_byc)) {

        /* Sanity check */
        assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_PE));

        const cs_real_t  *flux_array =
          input->values + cs_cdo_connect->c2e->idx[cm->c_id];
        for (short int e = 0; e < cm->n_ec; e++)
          fluxes[e] = flux_array[e];

      }
      else if (cs_test_flag(input->loc, cs_cdo_primal_cell)) {

        /* Retrieve the advection field:
           Switch to a cs_nvec3_t representation */
        cs_nvec3_t  adv_vect;
        cs_nvec3(input->values + 3*cm->c_id, &adv_vect);

        /* Sanity check */
        assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PEQ));

        /* Loop on cell edges */
        for (short int e = 0; e < cm->n_ec; e++)
          fluxes[e] = adv_vect.meas*cm->dface[e].meas*_dp3(adv_vect.unitv,
                                                           cm->dface[e].unitv);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at the cell center of cell %d.", adv->name, cm->c_id);
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)def->input;

      if (f->location_id == cs_mesh_location_get_id_by_name(N_("cells"))) {

        /* Retrieve the advection field:
           Switch to a cs_nvec3_t representation */
        cs_nvec3_t  adv_vect;
        cs_nvec3(f->val + 3*cm->c_id, &adv_vect);

        /* Sanity check */
        assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PEQ));

        /* Loop on cell edges */
        for (short int e = 0; e < cm->n_ec; e++)
          fluxes[e] = adv_vect.meas*cm->dface[e].meas*_dp3(adv_vect.unitv,
                                                           cm->dface[e].unitv);

      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Incompatible type of definition.");
    break;

  } // def_type

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         triangle defined by the two vertices of an edge and the barycenter
 *         of a face.
 *
 * \param[in]  adv       pointer to a cs_adv_field_t structure
 * \param[in]  cm        pointer to a cs_cell_mesh_t structure
 * \param[in]  tef_meas  area of the triangle tef
 * \param[in]  f         id of the face in the current cell
 * \param[in]  e         id of the edge in the current cell
 * \param[in]  v1        id of the first vertex in the current cell
 * \param[in]  v2        id of the second vertex in the current cell
 *
 * \return the value of the flux across tef
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_advection_field_get_flux_tef(const cs_adv_field_t        *adv,
                                const cs_cell_mesh_t        *cm,
                                const cs_real_t              tef_meas,
                                short int                    f,
                                short int                    e,
                                short int                    v1,
                                short int                    v2)
{
  cs_real_t  adv_flx = 0;

  if (adv == NULL)
    return adv_flx;

  /* Sanity check */
  assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ));

  cs_xdef_t  *def = adv->definition;

  const cs_real_t  *xv1 = cm->xv + 3*v1;
  const cs_real_t  *xv2 = cm->xv + 3*v2;
  const cs_quant_t  pfq = cm->face[f];

  /* Compute the flux accros the portion of primal face */
  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      /* Retrieve the advection field: Switch to a cs_nvec3_t representation */
      cs_real_3_t  cell_vector;
      cs_xdef_eval_cw_vector_by_val(cm, cs_time_step, def->input, cell_vector);
      cs_nvec3_t  adv_vect;
      cs_nvec3(cell_vector, &adv_vect);

      adv_flx = tef_meas * adv_vect.meas * _dp3(adv_vect.unitv, pfq.unitv);
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_real_t  eval[9];
      cs_quadrature_type_t  qtype = cs_xdef_get_quadrature(def);

      switch (qtype) {

      case CS_QUADRATURE_NONE:
      case CS_QUADRATURE_BARY:
      case CS_QUADRATURE_BARY_SUBDIV:
        {
          cs_real_3_t  xg;

          for (int k = 0; k < 3; k++)
            xg[k] = cs_math_onethird * (xv1[k] + xv2[k] + pfq.center[k]);

          /* Call the analytic function. result is stored in eval */
          cs_xdef_eval_cw_at_xyz_by_analytic(cm,
                                             1, (const cs_real_t *)xg,
                                             cs_time_step,
                                             def->input,
                                             eval);

          adv_flx = tef_meas * _dp3(eval, pfq.unitv);
        }
        break;

      case CS_QUADRATURE_HIGHER:
        {
          cs_real_t  w, add = 0.;
          cs_real_3_t  gpts[3];

          cs_quadrature_tria_3pts(xv1, xv2, pfq.center, tef_meas, gpts, &w);

          /* Call the analytic function. result is stored in eval for the three
             quadrature points */
          cs_xdef_eval_cw_at_xyz_by_analytic(cm,
                                             3, (const cs_real_t *)gpts,
                                             cs_time_step,
                                             def->input,
                                             eval);

          for (int p = 0; p < 3; p++) add += _dp3(eval + 3*p, pfq.unitv);
          adv_flx += add * w;

        }
        break;

      case CS_QUADRATURE_HIGHEST: // Not yet implemented
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of quadrature for computing the flux of %s"
                  " across an elementary triangle s(v,e,f).\n"
                  " This functionality is not implemented yet.", adv->name);
        break;
      }

    }
    break; // DEF_ANALYTIC_FUNCTION

    case CS_XDEF_BY_ARRAY:
      {
        cs_real_3_t  rec_field;
        cs_xdef_array_input_t  *input = (cs_xdef_array_input_t *)def->input;

        /* Test if flag has at least the pattern of the reference support */
        if (cs_test_flag(input->loc, cs_cdo_dual_face_byc)) {

          const cs_connect_index_t  *c2e = cs_cdo_connect->c2e;
          const cs_real_t  *cell_array = input->values + c2e->idx[cm->c_id];

          /* Compute the reconstruction of the flux in pec */
          cs_reco_dfbyc_in_pec(cm, e, cell_array, rec_field);

          /* The reconstruction yields a constant vector field */
          adv_flx = tef_meas * _dp3(pfq.unitv, rec_field);

        }
        else
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid support for evaluating the advection field %s"
                    " across tef.", adv->name);

      } // DEF_BY_ARRAY
      break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of definition for computing the flux of %s"
              " across the triangle tef.\n"
              " This functionality is not implemented yet.", adv->name);
    break;

  } // switch def_type

  return adv_flx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For each cs_adv_field_t structures, update the values of the related
 *         field(s)
 *
 * \param[in]      cur2prev    true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_update(bool   cur2prev)
{
  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];
    assert(adv != NULL);

    if (adv->vtx_field_id > -1) { /* Field stored at vertices */

      cs_field_t  *fld = cs_field_by_id(adv->vtx_field_id);

      /* Copy current field values to previous values */
      if (cur2prev)
        cs_field_current_to_previous(fld);

      /* Set the new values */
      cs_advection_field_at_vertices(adv, fld->val);

    }

    if (adv->cell_field_id > -1) { /* Field stored at cell centers */

      cs_field_t  *fld = cs_field_by_id(adv->cell_field_id);

      /* Copy current field values to previous values */
      if (cur2prev)
        cs_field_current_to_previous(fld);

      /* Set the new values */
      cs_advection_field_at_cells(adv, fld->val);

    }

  } // Loop on advection fields

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell
 *
 * \param[in]      adv        pointer to the advection field struct.
 * \param[in]      diff       pointer to the diffusion property struct.
 * \param[in, out] peclet     pointer to an array storing Peclet number
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_get_peclet(const cs_adv_field_t        *adv,
                        const cs_property_t         *diff,
                        cs_real_t                    peclet[])
{
  cs_real_t  ptymat[3][3];
  cs_real_3_t  ptydir;
  cs_nvec3_t  adv_c;

  const bool  pty_uniform = cs_property_is_uniform(diff);
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  assert(peclet != NULL); // Sanity check

  /* Get the value of the material property at the first cell center */
  if (pty_uniform)
    cs_property_get_cell_tensor(0, diff, false, ptymat);

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    /* Get the value of the material property at the cell center */
    if (!pty_uniform)
      cs_property_get_cell_tensor(c_id, diff, false, ptymat);

    cs_advection_field_get_cell_vector(c_id, adv, &adv_c);

    cs_real_t  hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);

    cs_math_33_3_product((const cs_real_t (*)[3])ptymat, adv_c.unitv, ptydir);

    peclet[c_id] = hc * adv_c.meas / _dp3(adv_c.unitv, ptydir);

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
