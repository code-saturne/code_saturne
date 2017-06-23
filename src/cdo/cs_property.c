/*============================================================================
 * Manage the definition/setting of properties
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_reco.h"
#include "cs_volume_zone.h"
#include "cs_xdef_eval.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_property.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_PROPERTY_DBG  0

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_pty[] =
  " Stop setting an empty cs_property_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

static int  _n_properties = 0;
static int  _n_max_properties = 0;
static cs_property_t  **_properties = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the settings are valid
 *
 * \param[in]  tens      values of the tensor
 */
/*----------------------------------------------------------------------------*/

static inline bool
_is_tensor_symmetric(const cs_real_3_t      *tens)
{
  if ((tens[0][1] - tens[1][0]) > cs_math_zero_threshold ||
      (tens[0][2] - tens[2][0]) > cs_math_zero_threshold ||
      (tens[1][2] - tens[2][1]) > cs_math_zero_threshold)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a tensor
 *
 * \param[in]  tensor      values of the tensor
 */
/*----------------------------------------------------------------------------*/

static inline void
_print_tensor(const cs_real_3_t   *tensor)
{
  cs_log_printf(CS_LOG_DEFAULT,
                "\n  Tensor property: | % 8.4e  % 8.4e  % 8.4e |\n"
                "                     | % 8.4e  % 8.4e  % 8.4e |\n"
                "                     | % 8.4e  % 8.4e  % 8.4e |\n",
                tensor[0][0], tensor[0][1], tensor[0][2],
                tensor[1][0], tensor[1][1], tensor[1][2],
                tensor[2][0], tensor[2][1], tensor[2][2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the settings are valid and then invert a tensor
 *
 * \param[in, out]  tens           values of the tensor
 * \param[in]       type           type of property to deal with
 */
/*----------------------------------------------------------------------------*/

static void
_invert_tensor(cs_real_3_t          *tens,
               cs_property_type_t    type)
{
#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 0 /* Sanity check */
  bool  is_ok =true;
  for (int k = 0; k < 3; k++)
    if (fabs(tensor[k][k]) < cs_math_zero_threshold)
      is_ok = false;

  if (is_ok)
    _is_tensor_symmetric(tens);

  if (!is_ok) {
    _print_tensor((const cs_real_t (*)[3])tens);
    bft_error(__FILE__, __LINE__, 0,
              " %s: A problem has been detected during the definition of the"
              " property %s in the cell %d.\n", __func__, pty->name, c_id);
  }
#endif

  if (type == CS_PROPERTY_ISO || type == CS_PROPERTY_ORTHO)
    for (int k = 0; k < 3; k++)
      tens[k][k] /= 1.0;

  else { /* anisotropic */

    cs_real_33_t  invmat;

    cs_math_33_inv_cramer((const cs_real_3_t (*))tens, invmat);
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        tens[k][l] = invmat[k][l];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new definition to a cs_property_t structure
 *         Sanity checks on the settings related to this definition.

 * \param[in, out]  pty       pointer to a cs_property_t structure
 *
 * \return the new definition id
 */
/*----------------------------------------------------------------------------*/

inline static int
_add_new_def(cs_property_t     *pty)
{
  int  new_id = pty->n_definitions;

  pty->n_definitions += 1;
  BFT_REALLOC(pty->defs, pty->n_definitions, cs_xdef_t *);
  BFT_REALLOC(pty->get_eval_at_cell, pty->n_definitions,
              cs_xdef_eval_t *);
  BFT_REALLOC(pty->get_eval_at_cell_cw, pty->n_definitions,
              cs_xdef_eval_cw_t *);

  return new_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name          name of the property
 * \param[in]  key_type      keyname of the type of property
 * \param[in]  n_subdomains  piecewise definition on n_subdomains
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_property_t *
_create_property(const char           *name,
                 int                   id,
                 cs_property_type_t    type)
{
  cs_property_t  *pty = NULL;

  BFT_MALLOC(pty, 1, cs_property_t);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, name, len);

  pty->id = id;
  pty->type = type;
  pty->state_flag = 0;

  pty->n_definitions = 0;
  pty->defs = NULL;
  pty->def_ids = NULL;

  pty->get_eval_at_cell = NULL;
  pty->get_eval_at_cell_cw = NULL;

  return pty;
}

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
cs_property_set_shared_pointers(const cs_cdo_quantities_t    *quant,
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
 * \brief  Retrieve the number of properties
 *
 * \return the number of properties
 */
/*----------------------------------------------------------------------------*/

int
cs_property_get_n_properties(void)
{
  return _n_properties;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name          name of the property
 * \param[in]  type          type of property
 * \param[in]  n_subdomains  piecewise definition on n_subdomains
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add(const char            *name,
                cs_property_type_t     type)
{
  cs_property_t  *pty = cs_property_by_name(name);

  if (pty != NULL) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" %s: An existing property has already the name %s.\n"
                    " Stop adding this property.\n"), __func__, name);
    return  pty;
  }

  int  pty_id = _n_properties;

  if (pty_id == 0) {
    _n_max_properties = 3;
    BFT_MALLOC(_properties, _n_max_properties, cs_property_t *);
  }

  _n_properties += 1;

  if (_n_properties > _n_max_properties) {
    _n_max_properties *= 2;
    BFT_REALLOC(_properties, _n_max_properties, cs_property_t *);
  }

  _properties[pty_id] = _create_property(name, pty_id, type);

  return _properties[pty_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the related property definition from its name
 *
 * \param[in]  name    name of the property to find
 *
 * \return NULL if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_name(const char   *name)
{
  if (_n_properties < 0)
    return NULL;
  assert(name != NULL);

  for (int i = 0; i < _n_properties; i++) {
    cs_property_t  *pty = _properties[i];
    assert(pty->name != NULL);
    if (strcmp(pty->name, name) == 0)
      return pty;
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the related property definition from its id
 *
 * \param[in]  id      id of the property to find
 *
 * \return NULL if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_id(int         id)
{
  if (_n_properties < 0)
    return NULL;
  if (id < 0 || id >= _n_max_properties)
    return NULL;

  return  _properties[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all cs_property_t structures and the array storing all the
 *         structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_destroy_all(void)
{
  if (_n_properties == 0)
    return;

  for (int i = 0; i < _n_properties; i++) {

    cs_property_t  *pty = _properties[i];

    if (pty == NULL)
      bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

    BFT_FREE(pty->name);
    BFT_FREE(pty->def_ids);

    for (int j = 0; j < pty->n_definitions; j++)
      pty->defs[j] = cs_xdef_free(pty->defs[j]);

    BFT_FREE(pty->defs);
    BFT_FREE(pty->get_eval_at_cell);
    BFT_FREE(pty->get_eval_at_cell_cw);

    BFT_FREE(pty);

  } // Loop on properties

  BFT_FREE(_properties);
  _n_properties = 0;
  _n_max_properties = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage of the definition of a property based on several
 *         definitions (i.e. definition by subdomains)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_finalize_setup(void)
{
  if (_n_properties == 0)
    return;

  for (int i = 0; i < _n_properties; i++) {

    cs_property_t  *pty = _properties[i];

    if (pty == NULL)
      bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

    if (pty->n_definitions > 1) { /* Initialization of def_ids */

      const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

      BFT_MALLOC(pty->def_ids, n_cells, short int);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < n_cells; j++)
        pty->def_ids[j] = -1; // Unset by default

      for (int id = 0; id < pty->n_definitions; id++) {

        const cs_xdef_t  *def = pty->defs[id];

        assert(def->z_id > 0);
        assert(def->support == CS_XDEF_SUPPORT_VOLUME);

        const cs_volume_zone_t  *z = cs_volume_zone_by_id(def->z_id);

        assert(z != NULL);

#       pragma omp parallel for if (z->n_cells > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < z->n_cells; j++)
          pty->def_ids[z->cell_ids[j]] = id;

      } // Loop on definitions

    }
    else { // n_definitions = 1

      if (pty->defs[0]->type == CS_XDEF_BY_VALUE)
        pty->state_flag |= CS_FLAG_STATE_UNIFORM;

    }

  } // Loop on properties

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the property is uniform, otherwise false
 *
 * \param[in]    pty    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_property_is_uniform(const cs_property_t   *pty)
{
  if (pty == NULL)
    return false;

  if (pty->state_flag & CS_FLAG_STATE_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of a property
 *
 * \param[in]    pty    pointer to a property
 *
 * \return  the name of the related property
 */
/*----------------------------------------------------------------------------*/

const char *
cs_property_get_name(const cs_property_t   *pty)
{
  if (pty == NULL)
    return NULL;

  return pty->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of a property
 *
 * \param[in]    pty    pointer to a property
 *
 * \return  the type of the related property
 */
/*----------------------------------------------------------------------------*/

cs_property_type_t
cs_property_get_type(const cs_property_t   *pty)
{
  if (pty == NULL)
    return CS_PROPERTY_N_TYPES;

  return pty->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an isotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty          pointer to a cs_property_t structure
 * \param[in]       zone_name    name of an already defined zone
 * \param[in]       val          value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_iso_by_value(cs_property_t    *pty,
                             const char       *zone_name,
                             double            val)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not isotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  const cs_volume_zone_t  *zone = cs_volume_zone_by_name(zone_name);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; // metadata
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        1, // dim
                                        zone->id,
                                        state_flag,
                                        meta_flag,
                                        &val);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_scalar_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_eval_cw_scalar_by_val;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an orthotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty         pointer to a cs_property_t structure
 * \param[in]       zone_name   name of an already defined zone
 * \param[in]       val         values to set (vector of size 3)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_ortho_by_value(cs_property_t    *pty,
                               const char       *zone_name,
                               double            val[])
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ORTHO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not orthotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  const cs_volume_zone_t  *zone = cs_volume_zone_by_name(zone_name);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; // metadata
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        3, // dim
                                        zone->id,
                                        state_flag,
                                        meta_flag,
                                        (void *)val);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_vector_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_eval_cw_vector_by_val;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty         pointer to a cs_property_t structure
 * \param[in]       zone_name   name of an already defined zone
 * \param[in]       tens        values to set (3x3 tensor)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_aniso_by_value(cs_property_t    *pty,
                               const char       *zone_name,
                               cs_real_t         tens[3][3])
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ANISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not anisotropic.\n"
              " Please check your settings.", pty->name);

  /* Check the symmetry */
  if (!_is_tensor_symmetric(tens))
    bft_error(__FILE__, __LINE__, 0,
              _(" The definition of the tensor related to the"
                " property %s is not symmetric.\n"
                " This case is not handled. Please check your settings.\n"),
              pty->name);

  int  new_id = _add_new_def(pty);
  const cs_volume_zone_t  *zone = cs_volume_zone_by_name(zone_name);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; // metadata
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        9, // dim
                                        zone->id,
                                        state_flag,
                                        meta_flag,
                                        tens);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_tensor_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_eval_cw_tensor_by_val;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function in
 *         a subdomain attached to the mesh location named ml_name
 *
 * \param[in, out]  pty         pointer to a cs_property_t structure
 * \param[in]       zone_name   name of an already defined zone
 * \param[in]       func        pointer to a cs_analytic_func_t function
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_analytic(cs_property_t        *pty,
                            const char           *zone_name,
                            cs_analytic_func_t   *func)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_def(pty);
  int  dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  const cs_volume_zone_t  *zone = cs_volume_zone_by_name(zone_name);
  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0; // metadata
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        dim,
                                        zone->id,
                                        state_flag,
                                        meta_flag,
                                        (void *)func);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_at_cells_by_analytic;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_eval_cw_cell_by_analytic;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to law depending on one
 *         scalar variable in a subdomain attached to the mesh location named
 *         ml_name
 *
 * \param[in, out] pty                  pointer to a cs_property_t structure
 * \param[in]      zone_name            name of an already defined zone
 * \param[in]      context              pointer to a structure (may be NULL)
 * \param[in]      get_eval_at_cell     pointer to a function (may be NULL)
 * \param[in]      get_eval_at_cell_cw  pointer to a function
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_func(cs_property_t         *pty,
                        const char            *zone_name,
                        void                  *context,
                        cs_xdef_eval_t        *get_eval_at_cell,
                        cs_xdef_eval_cw_t     *get_eval_at_cell_cw)
{
  int dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  int  def_id = _add_new_def(pty);
  const cs_volume_zone_t  *zone = cs_volume_zone_by_name(zone_name);
  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0; // metadata
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_FUNCTION,
                                        dim,
                                        zone->id,
                                        state_flag,
                                        meta_flag,
                                        context);

  pty->defs[def_id] = d;
  pty->get_eval_at_cell[def_id] = get_eval_at_cell;
  pty->get_eval_at_cell_cw[def_id] = get_eval_at_cell_cw;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_array(cs_property_t    *pty,
                         cs_flag_t         loc,
                         cs_real_t        *array,
                         cs_lnum_t        *index)
{
  int  id = _add_new_def(pty);
  assert(id == 0);
  /* zone->id = 0 since all the support is selected in this case */

  int dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  if (pty->n_definitions > 1)
    bft_error(__FILE__, __LINE__, 0,
              " When a definition by array is requested, the max. number"
              " of subdomains to consider should be equal to 1.\n"
              " Current value is %d for property %s.\n"
              " Please modify your settings.",
              pty->n_definitions, pty->name);

  cs_flag_t  state_flag = 0; // Will be updated during the creation
  cs_flag_t  meta_flag = 0;  // metadata
  cs_xdef_array_input_t  input = {.stride = dim,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        dim,
                                        0, // zone_id
                                        state_flag,
                                        meta_flag,
                                        &input);

  /* Set pointers */
  pty->defs[id] = d;

  if (dim == 1)
    pty->get_eval_at_cell[id] = cs_xdef_eval_scalar_at_cells_by_array;
  else
    pty->get_eval_at_cell[id] = cs_xdef_eval_nd_at_cells_by_array;
  pty->get_eval_at_cell_cw[id] = cs_xdef_eval_cw_cell_by_array;

  if (cs_test_flag(loc, cs_cdo_primal_cell)   == false &&
      cs_test_flag(loc, cs_cdo_primal_vtx)    == false &&
      cs_test_flag(loc, cs_cdo_dual_face_byc) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: case not available.\n", __func__);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_field(cs_property_t    *pty,
                         cs_field_t       *field)
{
  int  id = _add_new_def(pty);

  int dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  /* Sanity checks */
  assert(dim == field->dim);
  assert(id == 0);
  /* zone->id = 0 since all the support is selected in this case */

  const cs_volume_zone_t  *z = cs_volume_zone_by_id(0);
  if (field->location_id != z->location_id)
    bft_error(__FILE__, __LINE__, 0,
              " Property defined by field requests that the field location"
              " is supported by cells\n"
              " Property %s\n", pty->name);
  if (pty->n_definitions > 1)
    bft_error(__FILE__, __LINE__, 0,
              " When a definition by array is requested, the max. number"
              " of subdomains to consider should be equal to 1.\n"
              " Current value is %d for property %s.\n"
              " Please modify your settings.",
              pty->n_definitions, pty->name);

  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; // metadata

  pty->defs[id] = cs_xdef_volume_create(CS_XDEF_BY_FIELD,
                                        dim,
                                        0, // zone_id
                                        state_flag,
                                        meta_flag,
                                        field);

  pty->get_eval_at_cell[id] = cs_xdef_eval_cell_by_field;
  pty->get_eval_at_cell_cw[id] = cs_xdef_eval_cw_cell_by_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the value of the property at each cell. Store the
 *         evaluation in the given array.
 *
 * \param[in]       pty      pointer to a cs_property_t structure
 * \param[in, out]  array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_cells(const cs_property_t    *pty,
                          cs_real_t              *array)
{
  assert(pty != NULL && array != NULL);

  for (int i = 0; i < pty->n_definitions; i++) {

    cs_xdef_t  *def = pty->defs[i];
    const cs_volume_zone_t  *z = cs_volume_zone_by_id(def->z_id);

    pty->get_eval_at_cell[i](z->n_cells,
                             z->cell_ids,
                             false, // without compact output
                             cs_glob_mesh,
                             cs_cdo_connect,
                             cs_cdo_quant,
                             cs_time_step,
                             def->input,
                             array);

  } // Loop on definitions

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *
 * \param[in]      c_id           id of the current cell
 * \param[in]      pty            pointer to a cs_property_t structure
 * \param[in]      do_inversion   true or false
 * \param[in, out] tensor         3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_cell_tensor(cs_lnum_t               c_id,
                            const cs_property_t    *pty,
                            bool                    do_inversion,
                            cs_real_3_t            *tensor)
{
  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  int  def_id = 0;
  if (pty->n_definitions > 1)
    def_id = pty->def_ids[c_id];

  assert(pty->get_eval_at_cell[def_id] != NULL);

  cs_xdef_t  *def = pty->defs[def_id];

  switch (pty->type) {

  case CS_PROPERTY_ISO:
    {
      double  eval;

      pty->get_eval_at_cell[def_id](1,
                                    &c_id,
                                    true,  // compact output
                                    cs_glob_mesh,
                                    cs_cdo_connect,
                                    cs_cdo_quant,
                                    cs_time_step,
                                    def->input,
                                    &eval);

      tensor[0][0] = tensor[1][1] = tensor[2][2] = eval;
    }
    break;

  case CS_PROPERTY_ORTHO:
    {
      double  eval[3];

      pty->get_eval_at_cell[def_id](1,
                                    &c_id,
                                    true,  // compact output
                                    cs_glob_mesh,
                                    cs_cdo_connect,
                                    cs_cdo_quant,
                                    cs_time_step,
                                    def->input,
                                    eval);

      for (int k = 0; k < 3; k++)
        tensor[k][k] = eval[k];
    }
    break;

  case CS_PROPERTY_ANISO:
    pty->get_eval_at_cell[def_id](1,
                                  &c_id,
                                  true,  // compact output
                                  cs_glob_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  cs_time_step,
                                  def->input,
                                  (cs_real_t *)tensor);
    break;

  default:
    assert(0);
    break;
  }

  if (do_inversion)
    _invert_tensor(tensor, pty->type);

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT, "\n pty: %s, cell_id: %d\n", pty->name, c_id);
  _print_tensor(tensor);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *
 * \param[in]   c_id           id of the current cell
 * \param[in]   pty            pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_get_cell_value(cs_lnum_t              c_id,
                           const cs_property_t   *pty)
{
  cs_real_t  result = 0;
  int  def_id = 0;

  if (pty == NULL)
    return result;

  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  if (pty->n_definitions > 1)
    def_id = pty->def_ids[c_id];

  assert(pty->get_eval_at_cell[def_id] != NULL);

  cs_xdef_t  *def = pty->defs[def_id];

  pty->get_eval_at_cell[def_id](1,
                                &c_id,
                                true, // compact output
                                cs_glob_mesh,
                                cs_cdo_connect,
                                cs_cdo_quant,
                                cs_time_step,
                                def->input,
                                &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      pty            pointer to a cs_property_t structure
 * \param[in]      do_inversion   true or false
 * \param[in, out] tensor         3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_tensor_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_property_t    *pty,
                           bool                    do_inversion,
                           cs_real_3_t            *tensor)
{
  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  int  def_id = 0;
  if (pty->n_definitions > 1)
    def_id = pty->def_ids[cm->c_id];

  assert(pty->get_eval_at_cell_cw[def_id] != NULL);

  cs_xdef_t  *def = pty->defs[def_id];

  switch (pty->type) {

  case CS_PROPERTY_ISO:
    {
      double  eval;

      pty->get_eval_at_cell_cw[def_id](cm, cs_time_step, def->input, &eval);

      tensor[0][0] = tensor[1][1] = tensor[2][2] = eval;
    }
    break;

  case CS_PROPERTY_ORTHO:
    {
      double  eval[3];

      pty->get_eval_at_cell_cw[def_id](cm, cs_time_step, def->input, eval);

      for (int k = 0; k < 3; k++)
        tensor[k][k] = eval[k];
    }
    break;

  case CS_PROPERTY_ANISO:
    pty->get_eval_at_cell_cw[def_id](cm, cs_time_step, def->input,
                                     (cs_real_t *)tensor);
    break;

  default:
    assert(0);
    break;
  }

  if (do_inversion)
    _invert_tensor(tensor, pty->type);

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT, "\n pty: %s, cell_id: %d\n",
                pty->name, cm->c_id);
  _print_tensor(tensor);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]   cm        pointer to a cs_cell_mesh_t structure
 * \param[in]   pty       pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_value_in_cell(const cs_cell_mesh_t   *cm,
                          const cs_property_t    *pty)
{
  cs_real_t  result = 0;

  if (pty == NULL)
    return result;

  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  int  def_id = 0;
  if (pty->n_definitions > 1)
    def_id = pty->def_ids[cm->c_id];
  cs_xdef_t  *def = pty->defs[def_id];

  assert(pty->get_eval_at_cell_cw[def_id] != NULL);
  pty->get_eval_at_cell_cw[def_id](cm, cs_time_step, def->input, &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Fourier number in each cell
 *
 * \param[in]      pty        pointer to the diffusive property struct.
 * \param[in]      dt         value of the current time step
 * \param[in, out] fourier    pointer to an array storing Fourier numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_fourier(const cs_property_t     *pty,
                        double                   dt,
                        cs_real_t                fourier[])
{
  assert(fourier != NULL); // Sanity check
  assert(dt > 0.);

  const bool  pty_uniform = cs_property_is_uniform(pty);
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  if (pty->type == CS_PROPERTY_ISO) {

    cs_real_t  ptyval = 0.;
    if (pty_uniform)
      ptyval = cs_property_get_cell_value(0, pty);

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      if (!pty_uniform)
        ptyval = cs_property_get_cell_value(c_id, pty);

      const cs_real_t  hc2 = pow(cdoq->cell_vol[c_id], 2*cs_math_onethird);

      fourier[c_id] = dt * ptyval / hc2;

    } // Loop on cells

  }
  else { // Property is orthotropic or anisotropic

    cs_real_t  eig_max, eig_ratio;
    cs_real_t  ptymat[3][3];

    /* Get the value of the material property at the first cell center */
    if (pty_uniform) {
      cs_property_get_cell_tensor(0, pty, false, ptymat);
      cs_math_33_eigen((const cs_real_t (*)[3])ptymat, &eig_ratio, &eig_max);
    }

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      /* Get the value of the material property at the cell center */
      if (!pty_uniform) {
        cs_property_get_cell_tensor(c_id, pty, false, ptymat);
        cs_math_33_eigen((const cs_real_t (*)[3])ptymat, &eig_ratio, &eig_max);
      }

      const cs_real_t  hc2 = pow(cdoq->cell_vol[c_id], 2*cs_math_onethird);

      fourier[c_id] = dt * eig_max / hc2;

    } // Loop on cells

  } // Type of property

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a summary of the settings for all defined cs_property_t
 *         structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_log_setup(void)
{
  if (_n_properties == 0)
    return;
  assert(_properties != NULL);

  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the definition of properties\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_properties             %d\n",
                _n_properties);

  for (int i = 0; i < _n_properties; i++) {

    bool  is_uniform = false, is_steady = true;
    const cs_property_t  *pty = _properties[i];

    if (pty->state_flag & CS_FLAG_STATE_UNIFORM)  is_uniform = true;
    if (pty->state_flag & CS_FLAG_STATE_STEADY) is_steady = true;

    cs_log_printf(CS_LOG_SETUP, "\n <pty> %s uniform [%s], steady [%s], ",
                  pty->name,
                  cs_base_strtf(is_uniform), cs_base_strtf(is_steady));

    switch(pty->type) {
    case CS_PROPERTY_ISO:
      cs_log_printf(CS_LOG_SETUP, "type: isotropic\n");
      break;
    case CS_PROPERTY_ORTHO:
      cs_log_printf(CS_LOG_SETUP, "type: orthotropic\n");
      break;
    case CS_PROPERTY_ANISO:
      cs_log_printf(CS_LOG_SETUP, "type: anisotropic\n");
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of property."));
      break;
    }

    cs_log_printf(CS_LOG_SETUP, "       %s> n_subdomains    %d\n",
                  pty->name, pty->n_definitions);

    for (int j = 0; j < pty->n_definitions; j++)
      cs_xdef_log(pty->defs[j]);

    cs_log_printf(CS_LOG_SETUP, " </pty>");
  } // Loop on properties

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
