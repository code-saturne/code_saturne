/*============================================================================
 * Manage the definition/setting of properties
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_param_cdo.h"
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

static int  _n_properties = 0;
static int  _n_max_properties = 0;
static cs_property_t  **_properties = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
              cs_xdef_cw_eval_t *);

  return new_id;
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

  if (type & CS_PROPERTY_ISO || type & CS_PROPERTY_ORTHO)
    for (int k = 0; k < 3; k++)
      tens[k][k] = 1.0/tens[k][k];

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
 * \brief  Compute the value of a property at the cell center
 *
 * \param[in]   c_id     id of the current cell
 * \param[in]   t_eval   physical time at which one evaluates the term
 * \param[in]   pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_cell_value(cs_lnum_t              c_id,
                cs_real_t              t_eval,
                const cs_property_t   *pty)
{
  int  def_id = 0;
  if (pty->n_definitions > 1) {
    assert(pty->def_ids != NULL);
    def_id = pty->def_ids[c_id];
    assert(def_id > -1);
  }

  assert(pty->get_eval_at_cell[def_id] != NULL);

  cs_xdef_t  *def = pty->defs[def_id];
  cs_real_t  result = 0;

  pty->get_eval_at_cell[def_id](1, &c_id, true, /* dense output */
                                cs_glob_mesh,
                                cs_cdo_connect,
                                cs_cdo_quant,
                                t_eval,
                                def->context,
                                &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]  cm        pointer to a cs_cell_mesh_t structure
 * \param[in]  pty       pointer to a cs_property_t structure
 * \param[in]  t_eval    physical time at which one evaluates the term
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_value_in_cell(const cs_cell_mesh_t   *cm,
               const cs_property_t    *pty,
               cs_real_t               t_eval)
{
  cs_real_t  result = 0;
  int  def_id = 0;
  if (pty->n_definitions > 1) {
    def_id = pty->def_ids[cm->c_id];
    assert(def_id > -1);
  }
  cs_xdef_t  *def = pty->defs[def_id];

  assert(pty->get_eval_at_cell_cw[def_id] != NULL);
  pty->get_eval_at_cell_cw[def_id](cm, t_eval, def->context, &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached to a property at the cell
 *         center
 *
 * \param[in]      c_id        id of the current cell
 * \param[in]      t_eval      physical time at which one evaluates the term
 * \param[in]      pty         pointer to a cs_property_t structure
 * \param[in, out] tensor      3x3 matrix
 */
/*----------------------------------------------------------------------------*/

static void
_get_cell_tensor(cs_lnum_t               c_id,
                 cs_real_t               t_eval,
                 const cs_property_t    *pty,
                 cs_real_t               tensor[3][3])
{
  int  def_id = 0;
  if (pty->n_definitions > 1) {
    def_id = pty->def_ids[c_id];
    assert(def_id > -1);
  }
  assert(pty->get_eval_at_cell[def_id] != NULL);

  cs_xdef_t  *def = pty->defs[def_id];

  if (pty->type & CS_PROPERTY_ISO) {

    double  eval;
    pty->get_eval_at_cell[def_id](1, &c_id, true,  /* dense output */
                                  cs_glob_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  &eval);

    tensor[0][0] = tensor[1][1] = tensor[2][2] = eval;

  }
  else if (pty->type & CS_PROPERTY_ORTHO) {

    double  eval[3];
    pty->get_eval_at_cell[def_id](1, &c_id, true,  /* dense output */
                                  cs_glob_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  eval);

    for (int k = 0; k < 3; k++)
      tensor[k][k] = eval[k];

  }
  else {

    assert(pty->type & CS_PROPERTY_ANISO);
    pty->get_eval_at_cell[def_id](1, &c_id, true,  /* dense output */
                                  cs_glob_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  (cs_real_t *)tensor);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached to a property at the cell
 *         center. Case of a property as the product of two other properties.
 *
 * \param[in]      c_id          id of the current cell
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

static void
_get_cell_tensor_by_property_product(cs_lnum_t               c_id,
                                     cs_real_t               t_eval,
                                     const cs_property_t    *pty,
                                     cs_real_t               tensor[3][3])
{
  assert(pty->related_properties != NULL);
  const cs_property_t  *a = pty->related_properties[0];
  const cs_property_t  *b = pty->related_properties[1];

  cs_real_t  tensor_a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  cs_real_t  tensor_b[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  /* Evaluates each property */
  _get_cell_tensor(c_id, t_eval, a, tensor_a);
  _get_cell_tensor(c_id, t_eval, b, tensor_b);

  /* Compute the product */
  if (pty->type & CS_PROPERTY_ISO) {
    /*  a and b are isotropic */
    tensor[0][0] = tensor[1][1] = tensor[2][2] = tensor_a[0][0]*tensor_b[0][0];
  }
  else if (pty->type & CS_PROPERTY_ORTHO) {
    for (int k = 0; k < 3; k++)
      tensor[k][k] = tensor_a[k][k]*tensor_b[k][k];
  }
  else {
    assert(pty->type & CS_PROPERTY_ANISO);
    /* tensor has been initialized by the calling function */
    cs_math_33_product_add(tensor_a, tensor_b, tensor);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]      cm            pointer to a cs_cell_mesh_t structure
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

static void
_tensor_in_cell(const cs_cell_mesh_t   *cm,
                const cs_property_t    *pty,
                cs_real_t               t_eval,
                cs_real_t               tensor[3][3])
{
  int  def_id = 0;
  if (pty->n_definitions > 1) {
    def_id = pty->def_ids[cm->c_id];
    assert(def_id > -1);
  }
  assert(pty->get_eval_at_cell_cw[def_id] != NULL);

  cs_xdef_t  *def = pty->defs[def_id];

  if (pty->type & CS_PROPERTY_ISO) {

    double  eval;
    pty->get_eval_at_cell_cw[def_id](cm, t_eval, def->context, &eval);
    tensor[0][0] = tensor[1][1] = tensor[2][2] = eval;

  }
  else if (pty->type & CS_PROPERTY_ORTHO) {

    double  eval[3];
    pty->get_eval_at_cell_cw[def_id](cm, t_eval, def->context, eval);
    for (int k = 0; k < 3; k++)
      tensor[k][k] = eval[k];

  }
  else {

    assert(pty->type & CS_PROPERTY_ANISO);
    pty->get_eval_at_cell_cw[def_id](cm, t_eval, def->context,
                                     (cs_real_t *)tensor);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center.
 *         Version using a cs_cell_mesh_t structure and with a property
 *         defined as the product of two existing ones.
 *
 * \param[in]      cm            pointer to a cs_cell_mesh_t structure
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

static void
_tensor_in_cell_by_property_product(const cs_cell_mesh_t   *cm,
                                    const cs_property_t    *pty,
                                    cs_real_t               t_eval,
                                    cs_real_t               tensor[3][3])
{
  assert(pty->related_properties != NULL);
  const cs_property_t  *a = pty->related_properties[0];
  const cs_property_t  *b = pty->related_properties[1];

  cs_real_t  tensor_a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  cs_real_t  tensor_b[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  /* Evaluates each property */
  _tensor_in_cell(cm, a, t_eval, tensor_a);
  _tensor_in_cell(cm, b, t_eval, tensor_b);

  /* Compute the product */
  if (pty->type & CS_PROPERTY_ISO) {
    /*  a and b are isotropic */
    tensor[0][0] = tensor[1][1] = tensor[2][2] = tensor_a[0][0]*tensor_b[0][0];
  }
  else if (pty->type & CS_PROPERTY_ORTHO) {
    for (int k = 0; k < 3; k++)
      tensor[k][k] = tensor_a[k][k]*tensor_b[k][k];
  }
  else {
    assert(pty->type & CS_PROPERTY_ANISO);
    tensor[0][0] = tensor[1][1] = tensor[2][2] = 0;
    cs_math_33_product_add(tensor_a, tensor_b, tensor);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a property as the product of two existing properties
 *
 * \param[in, out]  pty      resulting property
 */
/*----------------------------------------------------------------------------*/

static void
_define_pty_by_product(cs_property_t          *pty)
{
  /* Only one definition is added in this case specifying that the definition
   * relies on other definitions to be defined. The exact way to specify values
   * is managed by the calling code (with a call to each sub-definition using
   * the standard algorithm)
   */

  int  id = _add_new_def(pty);
  assert(id == 0);

  int dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_SUB_DEFINITIONS,
                                        dim,
                                        0,     /* zone_id = all cells */
                                        state_flag,
                                        meta_flag,
                                        NULL); /* no input */

  /* Set pointers */
  pty->defs[id] = d;
  pty->get_eval_at_cell[id] = NULL;
  pty->get_eval_at_cell_cw[id] = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name       name of the property
 * \param[in]  id         id of the property to create
 * \param[in]  type       type of property
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_property_t *
_create_property(const char           *name,
                 int                   id,
                 cs_property_type_t    type)
{
  /* Check the sanity of type */
  if (type & CS_PROPERTY_ISO) {
    if (type & CS_PROPERTY_ANISO)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Detection of a wrong type for property %s\n"
                "Set to CS_PROPERTY_ISO and CS_PROPERTY_ANISO.",
                __func__, name);
    if (type & CS_PROPERTY_ORTHO)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Detection of a wrong type for property %s\n"
                "Set to CS_PROPERTY_ISO and CS_PROPERTY_ORTHO.",
                __func__, name);
  }
  else if (type & CS_PROPERTY_ORTHO) {
    if (type & CS_PROPERTY_ANISO)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Detection of a wrong type for property %s\n"
                "Set to CS_PROPERTY_ORTHO and CS_PROPERTY_ANISO.",
                __func__, name);
  }
  else {
    if ((type & CS_PROPERTY_ANISO) == 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: No type specified for property %s\n"
                " Set one among CS_PROPERTY_ISO, CS_PROPERTY_ORTHO or"
                " CS_PROPERTY_ANISO.", __func__, name);
  }

  cs_property_t  *pty = NULL;

  BFT_MALLOC(pty, 1, cs_property_t);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, name, len);

  pty->id = id;
  pty->type = type;
  pty->state_flag = 0;
  pty->process_flag = 0;

  pty->ref_value = 1.0;         /* default setting */

  pty->n_definitions = 0;
  pty->defs = NULL;
  pty->def_ids = NULL;

  pty->get_eval_at_cell = NULL;
  pty->get_eval_at_cell_cw = NULL;

  pty->n_related_properties = 0;
  pty->related_properties = NULL;

  return pty;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
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
 * \brief  Define a cs_property_t structure thanks to the product of two
 *         properties
 *         The type is infered from that of the related properties
 *         The value of the property is given as:
 *         value_ab = value_a * value_b
 *
 * \param[in]       name      name of the property
 * \param[in]       pty_a     pointer to a cs_property_t structure
 * \param[in]       pty_b     pointer to a cs_property_t structure
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add_as_product(const char             *name,
                           const cs_property_t    *pty_a,
                           const cs_property_t    *pty_b)
{
  if (pty_a == NULL || pty_b == NULL)
    return NULL;

  /* Determine the type of the new property */
  cs_property_type_t  type = CS_PROPERTY_BY_PRODUCT;

  /* pty_b  |pty_a -> iso | ortho | aniso
   * iso    |    iso      | ortho | aniso
   * ortho  |   ortho     | ortho | aniso
   * aniso  |   aniso     | aniso | aniso */

  if (pty_a->type & CS_PROPERTY_ISO) {
    if (pty_b->type & CS_PROPERTY_ISO)
      type |= CS_PROPERTY_ISO;
    else if (pty_b->type & CS_PROPERTY_ORTHO)
      type |= CS_PROPERTY_ORTHO;
    else if (pty_b->type & CS_PROPERTY_ANISO)
      type |= CS_PROPERTY_ANISO;
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of property.", __func__);
  }
  else if (pty_a->type & CS_PROPERTY_ANISO) {
    type |= CS_PROPERTY_ANISO;
  }
  else if (pty_a->type & CS_PROPERTY_ORTHO) {
    if (pty_b->type & CS_PROPERTY_ANISO)
      type |= CS_PROPERTY_ANISO;
    else
      type |= CS_PROPERTY_ORTHO;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of property.", __func__);

  cs_property_t  *pty_ab = cs_property_add(name, type);

  pty_ab->n_related_properties = 2;
  BFT_MALLOC(pty_ab->related_properties, 2, const cs_property_t *);

  pty_ab->related_properties[0] = pty_a;
  pty_ab->related_properties[1] = pty_b;

  return pty_ab;
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
 * \brief  Set optional parameters related to a cs_property_t structure
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       key       key related to a setting option
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_option(cs_property_t       *pty,
                       cs_property_key_t    key)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  switch(key) {

  case CS_PTYKEY_POST_FOURIER:
    pty->process_flag |= CS_PROPERTY_POST_FOURIER;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key not implemented for setting a property."));
    break;

  } /* Switch on keys */

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the reference value associated to a \ref cs_property_t structure
 *         This is a real number even whatever the type of property is.
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       refval   value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_reference_value(cs_property_t    *pty,
                                double            refval)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->ref_value = refval;
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

    if (pty->n_related_properties > 0)
      BFT_FREE(pty->related_properties);

    BFT_FREE(pty);

  } /* Loop on properties */

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

    if (pty->type & CS_PROPERTY_BY_PRODUCT)
      continue;

    if (pty->n_definitions > 1) { /* Initialization of def_ids */

      const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

      BFT_MALLOC(pty->def_ids, n_cells, short int);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < n_cells; j++)
        pty->def_ids[j] = -1; /* Unset by default */

      for (int id = 0; id < pty->n_definitions; id++) {

        const cs_xdef_t  *def = pty->defs[id];

        assert(def->z_id > 0);
        assert(def->support == CS_XDEF_SUPPORT_VOLUME);

        const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

        assert(z != NULL);

#       pragma omp parallel for if (z->n_elts > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < z->n_elts; j++)
          pty->def_ids[z->elt_ids[j]] = id;

      } /* Loop on definitions */

      /* Check if the property is defined everywhere */
      for (cs_lnum_t j = 0; j < n_cells; j++)
        if (pty->def_ids[j] == -1)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: cell %ld is unset for property %s\n",
                    __func__, (long)j, pty->name);

    }
    else if (pty->n_definitions == 0) {

      /* Default initialization */
      if (pty->type & CS_PROPERTY_ISO)
        cs_property_def_iso_by_value(pty, NULL, pty->ref_value);
      else if (pty->type & CS_PROPERTY_ORTHO) {
        cs_real_t  ref[3] =  {pty->ref_value, pty->ref_value, pty->ref_value};
        cs_property_def_ortho_by_value(pty, NULL, ref);
      }
      else if (pty->type & CS_PROPERTY_ANISO) {
        cs_real_t  ref[3][3] = { {pty->ref_value, 0, 0},
                                 {0, pty->ref_value, 0},
                                 {0, 0, pty->ref_value} };
        cs_property_def_aniso_by_value(pty, NULL, ref);
      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: Incompatible property type.",
                  __func__);

      cs_log_printf(CS_LOG_DEFAULT,
                    "\n Property \"%s\" will be defined using its reference"
                    " value.\n", pty->name);

    }

  } /* Loop on properties */

  for (int i = 0; i < _n_properties; i++) {

    cs_property_t  *pty = _properties[i];

    if (pty->type & CS_PROPERTY_BY_PRODUCT) {

      assert(pty->n_related_properties == 2);

      const cs_property_t  *pty_a = pty->related_properties[0];
      const cs_property_t  *pty_b = pty->related_properties[1];

      pty->ref_value = pty_a->ref_value * pty_b->ref_value;

      _define_pty_by_product(pty);

    } /* Only properties defined as a product */

  } /* Loop on properties */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_property_data_t structure. If property is NULL
 *         then one considers that this is a unitary property
 *
 * \param[in]      need_tensor  true if one needs a tensor-valued evaluation
 * \param[in]      need_eigen   true if one needs an evaluation of eigen values
 * \param[in]      property     pointer to the \ref cs_property_t structure
 * \param[in, out] data         structure to initialize (already allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_data_init(bool                     need_tensor,
                      bool                     need_eigen,
                      const cs_property_t     *property,
                      cs_property_data_t      *data)
{
  if (data == NULL)
    return;

  data->property = property;

  data->is_unity = false;
  data->is_iso = false;

  if (property == NULL) {
    data->is_iso = true;
    data->is_unity = true;
  }
  else {

    if (property->type & CS_PROPERTY_ISO) {
      data->is_iso = true;

      if (property->n_definitions == 1) {
        cs_xdef_t  *d = property->defs[0];
        if (d->type == CS_XDEF_BY_VALUE) {
          double  *dval = (double *)d->context;
          if (fabs(dval[0] - 1) < FLT_MIN)
            data->is_iso = true;
        }
      }
    }

  }

  const cs_real_t  ref_val = (property == NULL) ? 1. : property->ref_value;

  data->need_eigen = need_eigen;
  data->eigen_max = ref_val;
  data->eigen_ratio = 1.0;

  data->need_tensor = need_tensor;

  data->value = ref_val;
  data->tensor[0][0] = ref_val, data->tensor[0][1] = 0, data->tensor[0][2] = 0;
  data->tensor[1][0] = 0, data->tensor[1][1] = ref_val, data->tensor[1][2] = 0;
  data->tensor[2][0] = 0, data->tensor[2][1] = 0, data->tensor[2][2] = ref_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a single uniform and steady isotropic definition for the
 *         given cs_property_t structure.
 *         This is a specialized variant of \ref cs_property_def_iso_by_value
 *         since several assumptions are satisfied.
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       val      value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_constant_value(cs_property_t    *pty,
                               double            val)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not isotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);

  if (new_id > 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid setting: property %s is assumed to be constant.\n"
              " Several definitions have been added.\n"
              " Please check your settings.", __func__, pty->name);

  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        1,     /* dim */
                                        0,     /* all cells */
                                        state_flag,
                                        meta_flag,
                                        &val); /* context */

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_scalar_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_scalar_by_val;

  /* Set the state flag */
  pty->state_flag |= CS_FLAG_STATE_CELLWISE | CS_FLAG_STATE_STEADY;
  pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  /* Set automatically the reference value if all cells are selected */
  cs_property_set_reference_value(pty, val);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an isotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]       val      value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_iso_by_value(cs_property_t    *pty,
                             const char       *zname,
                             double            val)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not isotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_get_vol_zone_id(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        1,     /* dim */
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &val); /* context */

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_scalar_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_scalar_by_val;

  /* Set the state flag */
  pty->state_flag |= CS_FLAG_STATE_CELLWISE | CS_FLAG_STATE_STEADY;
  if (z_id == 0)
    pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  /* Set automatically the reference value if all cells are selected */
  if (z_id == 0)
    cs_property_set_reference_value(pty, val);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an orthotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]       val      values to set (vector of size 3)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_ortho_by_value(cs_property_t    *pty,
                               const char       *zname,
                               double            val[])
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ORTHO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not orthotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_get_vol_zone_id(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE |
        CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        3, /* dim */
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        (void *)val);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_vector_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_vector_by_val;

  /* Set the state flag */
  pty->state_flag |= CS_FLAG_STATE_CELLWISE | CS_FLAG_STATE_STEADY;
  if (z_id == 0)
    pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]       tens     values to set (3x3 tensor)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_aniso_by_value(cs_property_t    *pty,
                               const char       *zname,
                               cs_real_t         tens[3][3])
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ANISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not anisotropic.\n"
              " Please check your settings.", pty->name);

  /* Check the symmetry */
  if (!_is_tensor_symmetric((const cs_real_t (*)[3])tens))
    bft_error(__FILE__, __LINE__, 0,
              _(" The definition of the tensor related to the"
                " property %s is not symmetric.\n"
                " This case is not handled. Please check your settings.\n"),
              pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_get_vol_zone_id(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE |
        CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        9, /* dim */
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        tens);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_tensor_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_tensor_by_val;

  /* Set the state flag */
  pty->state_flag |= CS_FLAG_STATE_CELLWISE | CS_FLAG_STATE_STEADY;
  if (z_id == 0)
    pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function in
 *         a subdomain attached to the mesh location named ml_name
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]       func     pointer to a cs_analytic_func_t function
 * \param[in]       input    NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_time_func(cs_property_t      *pty,
                             const char         *zname,
                             cs_time_func_t     *func,
                             void               *input)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_def(pty);
  int  z_id = cs_get_vol_zone_id(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_time_func_context_t  tfc = { .func = func,
                                       .input = input,
                                       .free_input = NULL };

  /* Default initialization */
  pty->get_eval_at_cell[new_id] = NULL;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_by_time_func;

  int  dim = 0;
  if (pty->type & CS_PROPERTY_ISO) {
    dim = 1;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_scalar_at_cells_by_time_func;
  }
  else if (pty->type & CS_PROPERTY_ORTHO) {
    dim = 3;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_vector_at_cells_by_time_func;
  }
  else if (pty->type & CS_PROPERTY_ANISO) {
    dim = 9;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_tensor_at_cells_by_time_func;
  }
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Incompatible property type.",
              __func__);

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_TIME_FUNCTION,
                                        dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &tfc);

  pty->defs[new_id] = d;

  /* Set the state flag */
  pty->state_flag |= CS_FLAG_STATE_CELLWISE;
  if (z_id == 0)
    pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function in
 *         a subdomain attached to the mesh location named ml_name
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]       func     pointer to a cs_analytic_func_t function
 * \param[in]       input    NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_analytic(cs_property_t        *pty,
                            const char           *zname,
                            cs_analytic_func_t   *func,
                            void                 *input)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  cs_flag_t  state_flag = 0, meta_flag = 0; /* metadata */
  int  z_id = cs_get_vol_zone_id(zname);
  cs_xdef_analytic_context_t  ac = { .z_id = z_id,
                                     .func = func,
                                     .input = input,
                                     .free_input = NULL };

  int  dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &ac);

  int  new_id = _add_new_def(pty);
  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_at_cells_by_analytic;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_by_analytic;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to law depending on one
 *         scalar variable in a subdomain attached to the mesh location named
 *         ml_name
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]       context              pointer to a structure (may be NULL)
 * \param[in]       get_eval_at_cell     pointer to a function (may be NULL)
 * \param[in]       get_eval_at_cell_cw  pointer to a function
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_func(cs_property_t         *pty,
                        const char            *zname,
                        void                  *context,
                        cs_xdef_eval_t        *get_eval_at_cell,
                        cs_xdef_cw_eval_t     *get_eval_at_cell_cw)
{
  int  def_id = _add_new_def(pty);
  int  z_id = cs_get_vol_zone_id(zname);
  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0; /* metadata */

  int dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_FUNCTION,
                                        dim,
                                        z_id,
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
 * \param[in]       is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                            (true or false)
 * \param[in]       index     optional pointer to the array index
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_array(cs_property_t    *pty,
                         cs_flag_t         loc,
                         cs_real_t        *array,
                         bool              is_owner,
                         cs_lnum_t        *index)
{
  int  id = _add_new_def(pty);

  if (pty->n_definitions > 1)
    bft_error(__FILE__, __LINE__, 0,
              " When a definition by array is requested, the max. number"
              " of subdomains to consider should be equal to 1.\n"
              " Current value is %d for property %s.\n"
              " Please modify your settings.",
              pty->n_definitions, pty->name);

  int  dim = 1;
  if (pty->type == CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type == CS_PROPERTY_ANISO)
    dim = 9;

  cs_flag_t  state_flag = 0; /* Will be updated during the creation */
  cs_flag_t  meta_flag = 0;  /* metadata */

  /* z_id = 0 since all the support is selected in this case */
  cs_xdef_array_context_t  input = { .z_id = 0,
                                     .stride = dim,
                                     .loc = loc,
                                     .values = array,
                                     . is_owner = is_owner,
                                     .index = index };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        dim,
                                        0, /* zone_id */
                                        state_flag,
                                        meta_flag,
                                        &input);

  /* Set pointers */
  pty->defs[id] = d;

  if (dim == 1)
    pty->get_eval_at_cell[id] = cs_xdef_eval_scalar_at_cells_by_array;
  else
    pty->get_eval_at_cell[id] = cs_xdef_eval_nd_at_cells_by_array;
  pty->get_eval_at_cell_cw[id] = cs_xdef_cw_eval_by_array;

  if (cs_flag_test(loc, cs_flag_primal_cell)   == false &&
      cs_flag_test(loc, cs_flag_primal_vtx)    == false &&
      cs_flag_test(loc, cs_flag_dual_face_byc) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: case not available.\n", __func__);

  /* Set the state flag */
  if (cs_flag_test(loc, cs_flag_primal_cell))
    pty->state_flag |= CS_FLAG_STATE_CELLWISE;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a field structure
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
  /* z_id = 0 since all the support is selected in this case */

  const cs_zone_t  *z = cs_volume_zone_by_id(0);
  if (field->location_id != z->location_id)
    bft_error(__FILE__, __LINE__, 0,
              " Property defined by field requests that the field location"
              " is supported by cells\n"
              " Property %s\n", pty->name);
  if (pty->n_definitions > 1)
    bft_error(__FILE__, __LINE__, 0,
              " When a definition by field is requested, the max. number"
              " of subdomains to consider should be equal to 1.\n"
              " Current value is %d for property %s.\n"
              " Please modify your settings.",
              pty->n_definitions, pty->name);

  cs_flag_t  state_flag = CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; /* metadata */

  pty->defs[id] = cs_xdef_volume_create(CS_XDEF_BY_FIELD,
                                        dim,
                                        0, /* zone_id */
                                        state_flag,
                                        meta_flag,
                                        field);

  pty->get_eval_at_cell[id] = cs_xdef_eval_cell_by_field;
  pty->get_eval_at_cell_cw[id] = cs_xdef_cw_eval_by_field;

  /* Set the state flag */
  pty->state_flag |= CS_FLAG_STATE_CELLWISE;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the value of the property at each cell. Store the
 *         evaluation in the given array.
 *
 * \param[in]       t_eval   physical time at which one evaluates the term
 * \param[in]       pty      pointer to a cs_property_t structure
 * \param[in, out]  array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_cells(cs_real_t               t_eval,
                          const cs_property_t    *pty,
                          cs_real_t              *array)
{
  assert(pty != NULL && array != NULL);

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    assert(pty->related_properties != NULL);
    const cs_property_t  *a = pty->related_properties[0];
    const cs_property_t  *b = pty->related_properties[1];

    if (pty->type & CS_PROPERTY_ISO) {

      cs_real_t  *val_a = NULL;
      BFT_MALLOC(val_a, quant->n_cells, cs_real_t);
      memset(val_a, 0, quant->n_cells*sizeof(cs_real_t));

      /* 1. Evaluates the property A */
      for (int i = 0; i < a->n_definitions; i++) {

        cs_xdef_t  *def = a->defs[i];
        const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

        a->get_eval_at_cell[i](z->n_elts,
                               z->elt_ids,
                               false, /* without dense output */
                               cs_glob_mesh,
                               cs_cdo_connect,
                               quant,
                               t_eval,
                               def->context,
                               val_a);

      } /* Loop on definitions */

      /* 2. Evaluates the property B and operates the product */
      for (int i = 0; i < b->n_definitions; i++) {

        cs_xdef_t  *def = b->defs[i];
        const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

        b->get_eval_at_cell[i](z->n_elts,
                               z->elt_ids,
                               false, /* without dense output */
                               cs_glob_mesh,
                               cs_cdo_connect,
                               quant,
                               t_eval,
                               def->context,
                               array);

        for (cs_lnum_t j = 0; j < z->n_elts; j++)
          array[z->elt_ids[j]] *= val_a[z->elt_ids[j]];

      } /* Loop on definitions */

      BFT_FREE(val_a);

    }
    else {

      if (a->type & CS_PROPERTY_ISO) {

        cs_real_t  *val_a = NULL;
        BFT_MALLOC(val_a, quant->n_cells, cs_real_t);
        memset(val_a, 0, quant->n_cells*sizeof(cs_real_t));

        /* 1. Evaluates the property A */
        for (int i = 0; i < a->n_definitions; i++) {

          cs_xdef_t  *def = a->defs[i];
          const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

          a->get_eval_at_cell[i](z->n_elts,
                                 z->elt_ids,
                                 false, /* without dense output */
                                 cs_glob_mesh,
                                 cs_cdo_connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 val_a);

        } /* Loop on definitions */

        int  b_dim = (b->type & CS_PROPERTY_ORTHO) ? 3 : 9;

        /* 2. Evaluates the property B and operates the product */
        for (int i = 0; i < b->n_definitions; i++) {

          cs_xdef_t  *def = b->defs[i];
          const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

          b->get_eval_at_cell[i](z->n_elts,
                                 z->elt_ids,
                                 false, /* without dense output */
                                 cs_glob_mesh,
                                 cs_cdo_connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 array);

          for (cs_lnum_t j = 0; j < z->n_elts; j++) {
            const cs_real_t  acoef = val_a[z->elt_ids[j]];
            cs_real_t  *_a = array + b_dim*z->elt_ids[j];
            for (int k = 0; k < b_dim; k++)
              _a[k] *= acoef;
          }

        } /* Loop on definitions */

        BFT_FREE(val_a);

      }
      else if (b->type & CS_PROPERTY_ISO) {

        cs_real_t  *val_b = NULL;
        BFT_MALLOC(val_b, quant->n_cells, cs_real_t);
        memset(val_b, 0, quant->n_cells*sizeof(cs_real_t));

        /* 1. Evaluates the property B */
        for (int i = 0; i < b->n_definitions; i++) {

          cs_xdef_t  *def = b->defs[i];
          const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

          b->get_eval_at_cell[i](z->n_elts,
                                 z->elt_ids,
                                 false, /* without dense output */
                                 cs_glob_mesh,
                                 cs_cdo_connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 val_b);

        } /* Loop on definitions */

        int  a_dim = (a->type & CS_PROPERTY_ORTHO) ? 3 : 9;

        /* 2. Evaluates the property A and operates the product */
        for (int i = 0; i < a->n_definitions; i++) {

          cs_xdef_t  *def = a->defs[i];
          const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

          a->get_eval_at_cell[i](z->n_elts,
                                 z->elt_ids,
                                 false, /* without dense output */
                                 cs_glob_mesh,
                                 cs_cdo_connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 array);

          for (cs_lnum_t j = 0; j < z->n_elts; j++) {
            const cs_real_t  bcoef = val_b[z->elt_ids[j]];
            cs_real_t  *_a = array + a_dim*z->elt_ids[j];
            for (int k = 0; k < a_dim; k++)
              _a[k] *= bcoef;
          }

        } /* Loop on definitions */

        BFT_FREE(val_b);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Case not handled yet.\n", __func__);

    } /* Either a or b is an isotropic property */

  }
  else { /* Simple case: One has to evaluate the property */

    if ((pty->type & CS_PROPERTY_ISO) && cs_property_is_constant(pty)) {

#     pragma omp parallel for if (cs_cdo_connect->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cs_cdo_connect->n_cells; i++)
        array[i] = pty->ref_value;

    }
    else {

      for (int i = 0; i < pty->n_definitions; i++) {

        cs_xdef_t  *def = pty->defs[i];
        const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

        pty->get_eval_at_cell[i](z->n_elts,
                                 z->elt_ids,
                                 false, /* without dense output */
                                 cs_glob_mesh,
                                 cs_cdo_connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 array);

      } /* Loop on definitions */

    } /* Not isotropic and not constant */

  } /* Not defined as the product of two existing properties */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *
 * \param[in]      c_id          id of the current cell
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      do_inversion  true or false
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_cell_tensor(cs_lnum_t               c_id,
                            cs_real_t               t_eval,
                            const cs_property_t    *pty,
                            bool                    do_inversion,
                            cs_real_t               tensor[3][3])
{
  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  if (pty->type & CS_PROPERTY_BY_PRODUCT)
    _get_cell_tensor_by_property_product(c_id, t_eval, pty, tensor);
  else
    _get_cell_tensor(c_id, t_eval, pty, tensor);

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
 * \param[in]   c_id     id of the current cell
 * \param[in]   t_eval   physical time at which one evaluates the term
 * \param[in]   pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_get_cell_value(cs_lnum_t              c_id,
                           cs_real_t              t_eval,
                           const cs_property_t   *pty)
{
  cs_real_t  result = 0;

  if (pty == NULL)
    return result;

  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", __func__, pty->name);

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    assert(pty->related_properties != NULL);
    const cs_real_t  result_a = _get_cell_value(c_id, t_eval,
                                                pty->related_properties[0]);
    const cs_real_t  result_b = _get_cell_value(c_id, t_eval,
                                                pty->related_properties[1]);

    return result_a * result_b;

  }
  else {

    if (cs_property_is_constant(pty))
      return pty->ref_value;
    else
      return _get_cell_value(c_id, t_eval, pty);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached to a property at the cell
 *         center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]      cm            pointer to a cs_cell_mesh_t structure
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in]      do_inversion  true or false
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_tensor_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_property_t    *pty,
                           cs_real_t               t_eval,
                           bool                    do_inversion,
                           cs_real_t               tensor[3][3])
{
  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  if (pty->type & CS_PROPERTY_BY_PRODUCT)
    _tensor_in_cell_by_property_product(cm, pty, t_eval, tensor);
  else
    _tensor_in_cell(cm, pty, t_eval, tensor);

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
 * \param[in]  cm        pointer to a cs_cell_mesh_t structure
 * \param[in]  pty       pointer to a cs_property_t structure
 * \param[in]  t_eval    physical time at which one evaluates the term
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_value_in_cell(const cs_cell_mesh_t   *cm,
                          const cs_property_t    *pty,
                          cs_real_t               t_eval)
{
  cs_real_t  result = 0;

  if (pty == NULL)
    return result;

  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {
    assert(pty->related_properties != NULL);
    const cs_real_t  result_a = _value_in_cell(cm, pty->related_properties[0],
                                               t_eval);
    const cs_real_t  result_b = _value_in_cell(cm, pty->related_properties[1],
                                               t_eval);
    return result_a * result_b;
  }
  else {

    if (cs_property_is_constant(pty))
      return pty->ref_value;
    else
      return _value_in_cell(cm, pty, t_eval);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Fourier number in each cell
 *
 * \param[in]      pty       pointer to the diffusive property struct.
 * \param[in]      t_eval    physical time at which one evaluates the term
 * \param[in]      dt        value of the current time step
 * \param[in, out] fourier   pointer to an array storing Fourier numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_fourier(const cs_property_t    *pty,
                        cs_real_t               t_eval,
                        double                  dt,
                        cs_real_t               fourier[])
{
  assert(fourier != NULL); /* Sanity check */
  assert(dt > 0.);

  const bool  pty_uniform = cs_property_is_uniform(pty);
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  if (pty->type & CS_PROPERTY_ISO) {

    cs_real_t  ptyval = 0.;
    if (pty_uniform)
      ptyval = cs_property_get_cell_value(0, t_eval, pty);

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      const cs_real_t  hc = cbrt(cdoq->cell_vol[c_id]);
      if (!pty_uniform)
        ptyval = cs_property_get_cell_value(c_id, t_eval, pty);

      fourier[c_id] = dt * ptyval / (hc*hc);

    }

  }
  else { /* Property is orthotropic or anisotropic */

    cs_real_t  eig_max, eig_ratio;
    cs_real_t  ptymat[3][3];

    /* Get the value of the material property at the first cell center */
    if (pty_uniform) {
      cs_property_get_cell_tensor(0, t_eval, pty, false, ptymat);
      cs_math_33_eigen((const cs_real_t (*)[3])ptymat, &eig_ratio, &eig_max);
    }

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      const cs_real_t  hc = cbrt(cdoq->cell_vol[c_id]);

      /* Get the value of the material property at the cell center */
      if (!pty_uniform) {
        cs_property_get_cell_tensor(c_id, t_eval, pty, false, ptymat);
        cs_math_33_eigen((const cs_real_t (*)[3])ptymat, &eig_ratio, &eig_max);
      }

      fourier[c_id] = dt * eig_max / (hc*hc);

    }

  } /* Type of property */

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

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the definition of properties\n");
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h1);

  char  prefix[256];

  for (int i = 0; i < _n_properties; i++) {

    bool  is_uniform = false, is_steady = true;
    const cs_property_t  *pty = _properties[i];

    if (pty == NULL)
      continue;
    assert(strlen(pty->name) < 200); /* Check that prefix is large enough */

    if (pty->state_flag & CS_FLAG_STATE_UNIFORM)  is_uniform = true;
    if (pty->state_flag & CS_FLAG_STATE_STEADY) is_steady = true;

    cs_log_printf(CS_LOG_SETUP, "\n  * %s | Uniform %s Steady %s\n", pty->name,
                  cs_base_strtf(is_uniform), cs_base_strtf(is_steady));
    cs_log_printf(CS_LOG_SETUP, "  * %s | Reference value  % -8.4e\n",
                  pty->name, pty->ref_value);

    if (pty->type & CS_PROPERTY_ISO)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Type: isotropic", pty->name);
    else if (pty->type & CS_PROPERTY_ORTHO)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Type: orthotropic", pty->name);
    else if (pty->type & CS_PROPERTY_ANISO)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Type: anisotropic", pty->name);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: Invalid type of property."),
                __func__);

    if (pty->type & CS_PROPERTY_BY_PRODUCT)
      cs_log_printf(CS_LOG_SETUP, " | by product\n");
    else
      cs_log_printf(CS_LOG_SETUP, "\n");

    cs_log_printf(CS_LOG_SETUP, "  * %s | Number of definitions: %d\n\n",
                  pty->name, pty->n_definitions);

    for (int j = 0; j < pty->n_definitions; j++) {
      sprintf(prefix, "        Definition %3d", j);
      cs_xdef_log(prefix, pty->defs[j]);
    }

  } /* Loop on properties */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
