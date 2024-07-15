/*============================================================================
 * Manage the definition/setting of properties
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "bft_mem.h"

#include "cs_array.h"
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
static const cs_mesh_t  *cs_mesh;

static int  _n_properties = 0;
static int  _n_max_properties = 0;
static cs_property_t **_properties       = nullptr;

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

static inline int
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
 * \brief  Add a new definition at the boundary to a cs_property_t structure
 *         Sanity checks on the settings related to this definition.

 * \param[in, out]  pty       pointer to a cs_property_t structure
 *
 * \return the new definition id
 */
/*----------------------------------------------------------------------------*/

static inline int
_add_new_b_def(cs_property_t     *pty)
{
  int  new_id = pty->n_b_definitions;

  pty->n_b_definitions += 1;
  BFT_REALLOC(pty->b_defs, pty->n_b_definitions, cs_xdef_t *);

  return new_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the associated definition id for the current cell id
 *
 * \param[in]  c_id      current cell id
 * \param[in]  pty       pointer to a cs_property_t structure
 *
 * \return the definition id
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_cell_def_id(const cs_lnum_t c_id,
                 const cs_property_t *pty)
{
  if (pty->n_definitions > 1) {

    assert(pty->def_ids != nullptr);
    assert(pty->def_ids[c_id] > -1);
    return pty->def_ids[c_id];
  }
  else
    return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of a property at the cell center
 *
 * \param[in] c_id     id of the current cell
 * \param[in] t_eval   physical time at which one evaluates the term
 * \param[in] pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_cell_value(const cs_lnum_t      c_id,
                const cs_real_t      t_eval,
                const cs_property_t *pty)
{
  int def_id = _get_cell_def_id(c_id, pty);

  assert(pty->get_eval_at_cell[def_id] != nullptr);

  cs_xdef_t *def    = pty->defs[def_id];
  cs_real_t  result = 0;

  pty->get_eval_at_cell[def_id](1,
                                &c_id,
                                true, /* dense output */
                                cs_mesh,
                                cs_cdo_connect,
                                cs_cdo_quant,
                                t_eval,
                                def->context,
                                &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of an analytical property at the face center
 *
 * \param[in] f_id     id of the current face
 * \param[in] t_eval   physical time at which one evaluates the term
 * \param[in] def      pointer to a cs_xdef_t structure
 *
 * \return the value of the property for the given face
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_face_value_analytic(cs_lnum_t f_id, cs_real_t t_eval, const cs_xdef_t *def)
{
  assert(def != nullptr && def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_context_t *cx = (cs_xdef_analytic_context_t *)def->context;
  assert(cx != nullptr);

  /* Evaluate the function for this time at face barycenter */

  const cs_real_t *f_xc = cs_quant_get_face_center(f_id, cs_cdo_quant);

  cs_real_t result;

  cx->func(t_eval, 1, nullptr, f_xc, true, cx->input, &result);

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of a property at the face center
 *
 * \param[in] f_id     id of the current face
 * \param[in] t_eval   physical time at which one evaluates the term
 * \param[in] pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given face
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_face_value(cs_lnum_t f_id, cs_real_t t_eval, const cs_property_t *pty)
{

  const cs_adjacency_t *f2c = cs_cdo_connect->f2c;
  assert(f2c != nullptr);

  const cs_lnum_t n_cells = f2c->idx[f_id + 1] - f2c->idx[f_id];
  assert(n_cells == 1 || n_cells == 2);

  const cs_lnum_t *c_ids = f2c->ids + f2c->idx[f_id];

  /* Boundary face */
  if (n_cells == 1) {
    const int def_id = _get_cell_def_id(c_ids[0], pty);

    const cs_xdef_t *def = pty->defs[def_id];

    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) {
      return _get_face_value_analytic(f_id, t_eval, def);
    }
    return _get_cell_value(c_ids[0], t_eval, pty);
  }
  else {

    /* If same analytical property */
    if (_get_cell_def_id(c_ids[0], pty) == _get_cell_def_id(c_ids[1], pty)) {
      const cs_xdef_t *def = pty->defs[_get_cell_def_id(c_ids[0], pty)];

      if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) {
        return _get_face_value_analytic(f_id, t_eval, def);
      }
    }

    /* Average of shared cells scaled by cell volume */
    cs_real_t result = 0., vol = 0.;
    for (cs_lnum_t i_c = 0; i_c < n_cells; i_c++) {
      const cs_lnum_t c_id  = c_ids[i_c];
      const cs_real_t vol_c = cs_cdo_quant->cell_vol[c_id];

      result += vol_c * _get_cell_value(c_id, t_eval, pty);
      vol += vol_c;
    }

    return result / vol;
  }
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
  int         def_id = _get_cell_def_id(cm->c_id, pty);
  cs_xdef_t  *def = pty->defs[def_id];

  assert(pty->get_eval_at_cell_cw[def_id] != nullptr);
  pty->get_eval_at_cell_cw[def_id](cm, t_eval, def->context, &result);

  return result;
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

  else { /* anisotropic (sym. or not) */

    cs_real_33_t  invmat;

    cs_math_33_inv_cramer((const cs_real_3_t (*))tens, invmat);
    for (int ki = 0; ki < 3; ki++)
      for (int kj = 0; kj < 3; kj++)
        tens[ki][kj] = invmat[ki][kj];

  }
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
  int         def_id = _get_cell_def_id(c_id, pty);
  cs_xdef_t  *def = pty->defs[def_id];

  assert(pty->get_eval_at_cell[def_id] != nullptr);

  if (pty->type & CS_PROPERTY_ISO) {

    double  eval;
    pty->get_eval_at_cell[def_id](1, &c_id, true,  /* dense output */
                                  cs_mesh,
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
                                  cs_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  eval);

    for (int k = 0; k < 3; k++)
      tensor[k][k] = eval[k];

  }
  else if (pty->type & CS_PROPERTY_ANISO_SYM) {

    double  eval[6];
    pty->get_eval_at_cell[def_id](1, &c_id, true,  /* dense output */
                                  cs_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  eval);

    /* Diag. values */

    tensor[0][0] = eval[0];
    tensor[1][1] = eval[1];
    tensor[2][2] = eval[2];

    /* Extra-diag. values */

    tensor[0][1] = tensor[1][0] = eval[3];
    tensor[0][2] = tensor[2][0] = eval[4];
    tensor[1][2] = tensor[2][1] = eval[5];

  }
  else {

    assert(pty->type & CS_PROPERTY_ANISO);
    pty->get_eval_at_cell[def_id](1, &c_id, true,  /* dense output */
                                  cs_mesh,
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
  assert(pty->related_properties != nullptr);
  const cs_property_t  *a = pty->related_properties[0];
  const cs_property_t  *b = pty->related_properties[1];

  cs_real_t  tensor_a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  cs_real_t  tensor_b[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  /* Evaluate each property */

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
  int         def_id = _get_cell_def_id(cm->c_id, pty);
  cs_xdef_t  *def = pty->defs[def_id];

  assert(pty->get_eval_at_cell_cw[def_id] != nullptr);

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
  else if (pty->type & CS_PROPERTY_ANISO_SYM) {

    double  eval[6];
    pty->get_eval_at_cell_cw[def_id](cm, t_eval, def->context, eval);

    /* Diag. values */

    tensor[0][0] = eval[0];
    tensor[1][1] = eval[1];
    tensor[2][2] = eval[2];

    /* Extra-diag. values */

    tensor[0][1] = tensor[1][0] = eval[3];
    tensor[0][2] = tensor[2][0] = eval[4];
    tensor[1][2] = tensor[2][1] = eval[5];

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
  assert(pty->related_properties != nullptr);
  const cs_property_t  *a = pty->related_properties[0];
  const cs_property_t  *b = pty->related_properties[1];

  cs_real_t  tensor_a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  cs_real_t  tensor_b[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  /* Evaluate each property */

  _tensor_in_cell(cm, a, t_eval, tensor_a);
  _tensor_in_cell(cm, b, t_eval, tensor_b);

  /* Compute the product */

  if (pty->type & CS_PROPERTY_ISO) {
    /* a and b are isotropic */
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
  if (pty->type & CS_PROPERTY_ORTHO)
    dim = 3;
  else if (pty->type & CS_PROPERTY_ANISO_SYM)
    dim = 6;
  else if (pty->type & CS_PROPERTY_ANISO)
    dim = 9;

  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0;

  cs_xdef_t *d = cs_xdef_volume_create(CS_XDEF_BY_SUB_DEFINITIONS,
                                       dim,
                                       0, /* zone_id = all cells */
                                       state_flag,
                                       meta_flag,
                                       nullptr); /* no input */

  /* Set pointers */

  pty->defs[id] = d;
  pty->get_eval_at_cell[id]    = nullptr;
  pty->get_eval_at_cell_cw[id] = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill an array of values from a reference on a list of elements
 *
 * \param[in]      pty_dim   dimension of the reference value
 * \param[in]      n_elts    number of elements
 * \param[in]      elt_ids   list of element ids or nullptr
 * \param[in]      ref_val   reference values
 * \param[in, out] array     array storing the result of the evaluation(s)
 */
/*----------------------------------------------------------------------------*/

static void
_assign_ref_value(int                 pty_dim,
                  cs_lnum_t           n_elts,
                  const cs_lnum_t    *elt_ids,
                  const cs_real_t     ref_val[],
                  cs_real_t          *array)
{
  switch (pty_dim) {

  case 1: /* Isotropic */
    cs_array_real_set_scalar_on_subset(n_elts, elt_ids, ref_val[0], array);
    break;

  case 3: /* Orthotropic */
    cs_array_real_set_vector_on_subset(n_elts, elt_ids, ref_val, array);
    break;

  case 9: /* Anisotropic */
    {
      cs_real_t  tens[3][3] = {{ref_val[0], ref_val[1], ref_val[2]},
                               {ref_val[3], ref_val[4], ref_val[5]},
                               {ref_val[6], ref_val[7], ref_val[8]}};
      cs_array_real_set_tensor_on_subset(n_elts, elt_ids, tens, array);
    }
    break;

  default:
    /* Include the anisotropic with symmetric storage (pty_dim = 6) */
    cs_array_real_set_value_on_subset(n_elts, pty_dim, elt_ids, ref_val, array);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of a property at all boundary faces when the
 *        definition stems from the values at adjacent cells
 *
 * \param[in]      pty       pointer to a property structure
 * \param[in]      def_idx   pointer to the index on definitions
 * \param[in]      cell_ids  pointer to the related list of cell ids
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] eval      array storing the result of the evaluation(s)
 */
/*----------------------------------------------------------------------------*/

static void
_evaluate_property_at_boundary_from_cells(const cs_property_t    *pty,
                                          const cs_lnum_t        *def_idx,
                                          const cs_lnum_t        *cell_ids,
                                          double                  t_eval,
                                          cs_real_t              *eval)
{
  if (pty == nullptr)
    return;
  if (eval == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: Property \"%s\". Empty array.",
              __func__, pty->name);

  const cs_lnum_t  n_b_faces = cs_mesh->n_b_faces;
  const cs_lnum_t  *b_f2c = cs_mesh->b_face_cells;

  if (pty->n_definitions == 1) {

    cs_xdef_t  *def = pty->defs[0];

    pty->get_eval_at_cell[0](n_b_faces,
                             b_f2c,
                             true, /* dense output */
                             cs_mesh,
                             cs_cdo_connect,
                             cs_cdo_quant,
                             t_eval,
                             def->context,
                             eval);

  }
  else {

    const cs_lnum_t  *count = def_idx + pty->n_definitions + 1;

    for (int def_id = 0; def_id < pty->n_definitions; def_id++) {

      cs_xdef_t  *def = pty->defs[def_id];

      pty->get_eval_at_cell[def_id](count[def_id],
                                    cell_ids + def_idx[def_id],
                                    true, /* dense output */
                                    cs_mesh,
                                    cs_cdo_connect,
                                    cs_cdo_quant,
                                    t_eval,
                                    def->context,
                                    eval + def_idx[def_id]);

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate if needed and build an index definition by definition of the
 *        list of associated boundary cells. Useful to evaluate the value of a
 *        property at the boundary from the knowledge of the property value
 *        associated to a cell having boundary faces
 *        If def_idx is allocated then its size should be 2*n_definitions+1
 *        If cell_ids is allocated then its size should be n_b_faces
 *        If bf_ids is allocated then its size should be n_b_faces
 *
 * \param[in]  pty         pointer to a property structure
 * \param[in]  p_def_idx   double pointer to the index on definitions
 * \param[in]  p_cell_ids  double pointer to the related list of cell ids
 * \param[in]  p_bf_ids    double pointer to the list of related boundary faces
 */
/*----------------------------------------------------------------------------*/

static void
_build_def_idx_from_bdy_selection(const cs_property_t    *pty,
                                  cs_lnum_t             **p_def_idx,
                                  cs_lnum_t             **p_cell_ids,
                                  cs_lnum_t             **p_bf_ids)
{
  if (pty == nullptr)
    return;
  if (pty->n_definitions == 1)
    return; /* This can be done easier */

  assert(pty->def_ids != nullptr);

  const cs_lnum_t  n_b_faces = cs_mesh->n_b_faces;
  const cs_lnum_t  *b_f2c = cs_mesh->b_face_cells;

  cs_lnum_t  *def_idx = *p_def_idx;
  cs_lnum_t  *cell_ids = *p_cell_ids;
  cs_lnum_t  *bf_ids = *p_bf_ids;

  /* Allocations and initializations */

  if (def_idx == nullptr)
    BFT_MALLOC(def_idx, 2*pty->n_definitions + 1, cs_lnum_t);

  memset(def_idx, 0, sizeof(cs_lnum_t)*(2*pty->n_definitions + 1));

  cs_lnum_t  *count = def_idx + pty->n_definitions + 1;

  if (cell_ids == nullptr)
    BFT_MALLOC(cell_ids, n_b_faces, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    cell_ids[i] = -1;

  if (bf_ids == nullptr)
    BFT_MALLOC(bf_ids, n_b_faces, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    bf_ids[i] = -1;

  /* Build the index */

  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    def_idx[pty->def_ids[b_f2c[i]] + 1] += 1;

  for (int i = 0; i < pty->n_definitions; i++)
    def_idx[i+1] += def_idx[i];

  /* Define the list of cell ids for each definition */

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {

    const cs_lnum_t  c_id = b_f2c[i];
    const short int  def_id = pty->def_ids[c_id];
    const cs_lnum_t  shift = def_idx[def_id] + count[def_id];

    cell_ids[shift] = c_id;
    bf_ids[shift] = i;
    count[def_id] += 1;

  }

  /* Return pointers */

  *p_def_idx = def_idx;
  *p_cell_ids = cell_ids;
  *p_bf_ids = bf_ids;
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
  int n_types = 0;
  const int flags[] = {CS_PROPERTY_ISO,
                       CS_PROPERTY_ORTHO,
                       CS_PROPERTY_ANISO_SYM,
                       CS_PROPERTY_ANISO};

  for (int i = 0; i < 4; i++) {
    if (type & flags[i])
      n_types += 1;
  }

  if (n_types > 1) {

    const char *names[] = {"CS_PROPERTY_ISO",
                           "CS_PROPERTY_ORTHO",
                           "CS_PROPERTY_ANISO_SYM",
                           "CS_PROPERTY_ANISO"};
    int l = 0;
    char prop_list[256] = "";

    for (int i = 0; i < 4 && l > 0; i++) {
      if (type & flags[i]) {
        snprintf(prop_list+l, 256-l, "  %s\n", names[i]);
        prop_list[255] = '\0';
        l = strlen(prop_list);
      }
    }

  }
  else if (n_types < 1)
    if ((type & CS_PROPERTY_ANISO) == 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: No known type specified for property %s\n"
                " Set one among\n"
                "   CS_PROPERTY_ISO,\n"
                "   CS_PROPERTY_ORTHO,\n"
                "   CS_PROPERTY_ANISO_SYM,\n"
                "   CS_PROPERTY_ANISO.\n", __func__, name);

  cs_property_t *pty = nullptr;
  BFT_MALLOC(pty, 1, cs_property_t);

  /* Copy name */

  size_t  len = strlen(name);
  BFT_MALLOC(pty->name, len + 1, char);
  strncpy(pty->name, name, len + 1);

  pty->id = id;
  pty->type = type;
  pty->state_flag = 0;
  pty->process_flag = 0;

  pty->ref_value = 1.0;         /* default setting */
  pty->scaling_factor = 1.0;    /* default setting */

  pty->n_definitions = 0;
  pty->defs          = nullptr;
  pty->def_ids       = nullptr;

  pty->get_eval_at_cell    = nullptr;
  pty->get_eval_at_cell_cw = nullptr;

  pty->n_related_properties = 0;
  pty->related_properties   = nullptr;

  pty->n_b_definitions = 0;
  pty->b_defs          = nullptr;
  pty->b_def_ids       = nullptr;

  return pty;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set shared pointers to main domain members
 *
 * \param[in] mesh        mesh structure shared between FV and CDO
 * \param[in] quant       additional mesh quantities struct.
 * \param[in] connect     pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_property_init_sharing(const cs_mesh_t              *mesh,
                         const cs_cdo_quantities_t    *quant,
                         const cs_cdo_connect_t       *connect)
{
  cs_mesh = mesh;
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the number of properties
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
 * \brief Create and initialize a new property structure
 *
 * \param[in] name        name of the property
 * \param[in] type        type of property
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add(const char            *name,
                cs_property_type_t     type)
{
  cs_property_t  *pty = cs_property_by_name(name);

  if (pty != nullptr) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
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
 * \brief Create and initialize a new property structure with an evaluation
 *        which can be called on a sub-partition of a cell.
 *        This kind of property is not available for all numerical scheme.
 *        By default, only one evaluation is performed in each cell.
 *
 * \param[in] name          name of the property
 * \param[in] type          type of property
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_subcell_add(const char            *name,
                        cs_property_type_t     type)
{
  return cs_property_add(name, type | CS_PROPERTY_SUBCELL_DEFINITION);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_property_t structure thanks to the product of two
 *        properties
 *        The type is infered from that of the related properties
 *        The value of the property is given as value_ab = value_a * value_b
 *
 * \param[in] name      name of the property
 * \param[in] pty_a     pointer to a cs_property_t structure
 * \param[in] pty_b     pointer to a cs_property_t structure
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add_as_product(const char             *name,
                           const cs_property_t    *pty_a,
                           const cs_property_t    *pty_b)
{
  if (pty_a == nullptr || pty_b == nullptr)
    return nullptr;

  /* Determine the type for the new property */

  cs_property_type_t  type = CS_PROPERTY_BY_PRODUCT;

  /*              | pty_a->iso | pty_a->ortho | pty_a->aniso
   * pty_b->iso   |    iso     |   ortho      |    aniso
   * pty_b->ortho |   ortho    |   ortho      |    aniso
   * pty_b->aniso |   aniso    |   aniso      |    aniso
   */

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
 * \brief Find the related property definition from its name
 *
 * \param[in] name    name of the property to find
 *
 * \return nullptr if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_name(const char   *name)
{
  if (_n_properties < 0)
    return nullptr;
  assert(name != nullptr);

  for (int i = 0; i < _n_properties; i++) {
    cs_property_t  *pty = _properties[i];
    assert(pty->name != nullptr);
    if (strcmp(pty->name, name) == 0)
      return pty;
  }

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find the related property definition from its id
 *
 * \param[in] id      id of the property to find
 *
 * \return nullptr if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_id(int         id)
{
  if (_n_properties < 0)
    return nullptr;
  if (id < 0 || id >= _n_max_properties)
    return nullptr;

  return  _properties[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set optional parameters related to a cs_property_t structure
 *
 * \param[in, out] pty       pointer to a cs_property_t structure
 * \param[in]      key       key related to a setting option
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_option(cs_property_t       *pty,
                       cs_property_key_t    key)
{
  if (pty == nullptr)
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
 * \brief Set the reference value associated to a \ref cs_property_t structure
 *        This is a real number even whatever the type of property is.
 *
 * \param[in, out] pty      pointer to a cs_property_t structure
 * \param[in]      refval   value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_reference_value(cs_property_t    *pty,
                                double            refval)
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->ref_value = refval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling factor associated to a \ref cs_property_t structure
 *        This is a real number whatever the type of property is. If the
 *        property was not defined as CS_PROPERTY_SCALED, then this tag is
 *        added.
 *
 * \param[in, out] pty    pointer to a cs_property_t structure
 * \param[in]      val    value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_scaling_factor(cs_property_t    *pty,
                               double            val)
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->scaling_factor = val;
  pty->type |= CS_PROPERTY_SCALED;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling factor associated to a \ref cs_property_t structure
 *        This is a real number whatever the type of property is. If the
 *        property was not defined as CS_PROPERTY_SCALED, then this tag is
 *        added.
 *
 * \param[in, out] pty    pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_unscale(cs_property_t    *pty)
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->scaling_factor = 1.0;
  if (pty->type & CS_PROPERTY_SCALED)
    pty->type -= CS_PROPERTY_SCALED;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all cs_property_t structures and the array storing all the
 *        structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_destroy_all(void)
{
  if (_n_properties == 0)
    return;

  for (int i = 0; i < _n_properties; i++) {

    cs_property_t  *pty = _properties[i];

    if (pty == nullptr)
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

    if (pty->n_b_definitions > 0) {

      BFT_FREE(pty->b_defs);
      if (pty->n_b_definitions > 1)
        BFT_FREE(pty->b_def_ids);

    }

    BFT_FREE(pty);

  } /* Loop on properties */

  BFT_FREE(_properties);
  _n_properties = 0;
  _n_max_properties = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last stage of the definition of a property based on several
 *        definitions (i.e. definition by subdomains)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_finalize_setup(void)
{
  if (_n_properties == 0)
    return;

  for (int i = 0; i < _n_properties; i++) {

    cs_property_t  *pty = _properties[i];

    if (pty == nullptr)
      bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

    if (pty->type & CS_PROPERTY_BY_PRODUCT)
      continue; /* This is done after */

    /* Volume definitions */
    /* ------------------ */

    if (pty->n_definitions > 1) { /* Initialization of def_ids */

      const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

      BFT_MALLOC(pty->def_ids, n_cells, short int);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < n_cells; j++)
        pty->def_ids[j] = -1; /* Unset by default */

      for (int id = 0; id < pty->n_definitions; id++) {

        cs_xdef_t  *def = pty->defs[id];

        assert(def->z_id > 0);
        assert(def->support == CS_XDEF_SUPPORT_VOLUME);

        const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
        assert(z != nullptr);

#       pragma omp parallel for if (z->n_elts > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < z->n_elts; j++)
          pty->def_ids[z->elt_ids[j]] = id;

        /* If the definition is by array on a subset with an array allocated on
         * this subset, then one has to define an additional array for handling
         * the indirection between the full set and the subset.
         */

        if (def->type == CS_XDEF_BY_ARRAY) {

          cs_xdef_array_context_t *cx
            = static_cast<cs_xdef_array_context_t *>(def->context);

          if (!cx->full_length) {

            if (def->z_id != cx->z_id)
              bft_error(__FILE__, __LINE__, 0,
                        "%s: Issue found with the volume definition by array"
                        " for the property \"%s\"\n", __func__, pty->name);

            cs_xdef_array_build_full2subset(def);

          }

        } /* Definition by array */

      } /* Loop on definitions */

      /* Check if a definition id is associated to each cell */

      for (cs_lnum_t j = 0; j < n_cells; j++)
        if (pty->def_ids[j] == -1)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: cell %ld is unset for the property \"%s\"\n",
                    __func__, (long)j, pty->name);

    }
    else if (pty->n_definitions == 0) {

      /* Default definition based on the reference value */

      if (pty->type & CS_PROPERTY_ISO)
        cs_property_def_iso_by_value(pty, nullptr, pty->ref_value);
      else if (pty->type & CS_PROPERTY_ORTHO) {
        cs_real_t  ref[3] =  {pty->ref_value, pty->ref_value, pty->ref_value};
        cs_property_def_ortho_by_value(pty, nullptr, ref);
      }
      else if (pty->type & CS_PROPERTY_ANISO) {
        cs_real_t  ref[3][3] = { {pty->ref_value, 0, 0},
                                 {0, pty->ref_value, 0},
                                 {0, 0, pty->ref_value} };
        cs_property_def_aniso_by_value(pty, nullptr, ref);
      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: Incompatible property type.",
                  __func__);

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "\n The property \"%s\" will be defined using its reference"
                    " value.\n", pty->name);

    }

    /* Boundary definitions */
    /* -------------------- */

    if (pty->n_b_definitions > 1) { /* Initialization of b_def_ids */

      const cs_lnum_t  n_b_faces = cs_cdo_quant->n_b_faces;

      BFT_MALLOC(pty->b_def_ids, n_b_faces, short int);

#     pragma omp parallel for if (n_b_faces > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < n_b_faces; j++)
        pty->b_def_ids[j] = -1; /* Unset by default */

      for (int id = 0; id < pty->n_b_definitions; id++) {

        cs_xdef_t  *def = pty->b_defs[id];

        assert(def->z_id > 0);
        assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

        const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
        assert(z != nullptr);

#       pragma omp parallel for if (z->n_elts > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < z->n_elts; j++)
          pty->b_def_ids[z->elt_ids[j]] = id;

        /* If the definition is by array on a subset with an array allocated on
         * this subset, then one has to define an additional array for handling
         * the indirection between the full set and the subset.
         */

        if (def->type == CS_XDEF_BY_ARRAY) {

          cs_xdef_array_context_t *cx
            = static_cast<cs_xdef_array_context_t *>(def->context);

          if (!cx->full_length) {

            if (def->z_id != cx->z_id)
              bft_error(__FILE__, __LINE__, 0,
                        "%s: Issue found with the boundary definition by array"
                        " for the property \"%s\"\n", __func__, pty->name);

            cs_xdef_array_build_full2subset(def);

          }

        } /* Definition by array */

      } /* Loop on definitions */

      /* Check if a definition id is associated to each boundary face */

      for (cs_lnum_t j = 0; j < n_b_faces; j++)
        if (pty->b_def_ids[j] == -1)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Boundary face %ld is unset for the property \"%s\"\n",
                    __func__, (long)j, pty->name);

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
 * \brief Retrieve the array associated to the volume definition for the given
 *        property.
 *        Available only if there is one definition by array for the volume.
 *
 * \param[in] pty     pointer to the property structure
 *
 * \return a pointer to the array or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_property_get_array(const cs_property_t     *pty)
{
  if (pty == nullptr)
    return nullptr;
  if (pty->n_definitions > 1)
    return nullptr; /* May be too restrictive */

  const cs_xdef_t  *def = pty->defs[0];
  assert(def != nullptr);

  if (def->type == CS_XDEF_BY_ARRAY)
    return cs_xdef_array_get_values(def);
  else
    return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a \ref cs_property_data_t structure (not a pointer to this
 *         structure). If property is nullptr then one considers that this is a
 *         unitary property
 *
 * \param[in]   need_tensor  true if one needs a tensor-valued evaluation
 * \param[in]   need_eigen   true if one needs an evaluation of eigen values
 * \param[in]   property     pointer to the \ref cs_property_t structure
 *
 * \return an initialized structure
 */
/*----------------------------------------------------------------------------*/

cs_property_data_t
cs_property_data_define(bool                     need_tensor,
                        bool                     need_eigen,
                        const cs_property_t     *property)
{
  cs_property_data_t data;

  // Default initialization for the first set of members

  data.property = property;
  data.is_iso = false;
  data.is_unity = false;
  data.need_tensor = need_tensor;
  data.need_eigen = need_eigen;
  data.eigen_ratio = 1.0; // Computed only if needed

  if (property == nullptr)
    data.is_iso = true, data.is_unity = true;

  else {

    if (property->type & CS_PROPERTY_ISO) {
      data.is_iso = true;

      if (property->n_definitions == 1) {
        cs_xdef_t  *d = property->defs[0];
        if (d->type == CS_XDEF_BY_VALUE) {
          double  *dval = (double *)d->context;
          if (fabs(dval[0] - 1) < FLT_MIN)
            data.is_unity = true;
        }
      }
    }

  } // property != nullptr

  const cs_real_t ref_val = (property == nullptr) ? 1. : property->ref_value;

  // Second set of members

  data.eigen_max = ref_val;
  data.value = ref_val;
  data.tensor[0][0] = ref_val, data.tensor[0][1] = 0, data.tensor[0][2] = 0;
  data.tensor[1][0] = 0, data.tensor[1][1] = ref_val, data.tensor[1][2] = 0;
  data.tensor[2][0] = 0, data.tensor[2][1] = 0, data.tensor[2][2] = ref_val;

  return data;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_property_data_t structure. If property is
 * nullptr then one considers that this is a unitary property
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
  if (data == nullptr)
    return;

  data->property = property;
  data->is_unity = false;
  data->is_iso = false;

  if (property == nullptr) {
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
            data->is_unity = true;
        }
      }
    }

  }

  const cs_real_t ref_val = (property == nullptr) ? 1. : property->ref_value;

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
  if (pty == nullptr)
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
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
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
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not isotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
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
 * \brief Define the value at the boundary for the given (isotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]       val      value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_iso_by_value(cs_property_t    *pty,
                                      const char       *zname,
                                      double            val)
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not isotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          1,     /* dim */
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          &val); /* context */

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an orthotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
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
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ORTHO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not orthotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
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
 * \brief Define the value at the boundary for the given (orthotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]       vals     values to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_ortho_by_value(cs_property_t    *pty,
                                        const char       *zname,
                                        double            vals[])
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ORTHO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not orthotropic.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          3,     /* dim */
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          vals); /* context */

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
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
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ANISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not anisotropic.\n"
              " Please check your settings.", pty->name);

  /* Check the symmetry */

  if (!_is_tensor_symmetric((const cs_real_t (*)[3])tens))
    bft_error(__FILE__, __LINE__, 0,
              " %s: The definition of the tensor related to the property"
              " \"%s\" is not symmetric.\n"
              " This case is not handled. Please check your settings.\n",
              __func__, pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
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
 * \brief Define the value at the boundary for the given (anisotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]       tens     values to set given as a 3x3 tensor
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_aniso_by_value(cs_property_t    *pty,
                                        const char       *zname,
                                        double            tens[3][3])
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ANISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not orthotropic.\n"
              " Please check your settings.", pty->name);
  if (!_is_tensor_symmetric((const cs_real_t (*)[3])tens))
    bft_error(__FILE__, __LINE__, 0,
              " %s: The definition of the tensor related to the property"
              " \"%s\" is not symmetric.\n"
              " This case is not handled. Please check your settings.\n",
              __func__, pty->name);

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          9,      /* dim */
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          tens);  /* context */

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure thanks to a time
 *        function for all cells associated to the zone named zname.
 *        Optimized case with a symmetric storage.
 *
 * \param[in, out] pty       pointer to a cs_property_t structure
 * \param[in]      zname     name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]      symtens   symmetric tensor given as an array of 6 values
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_aniso_sym_by_value(cs_property_t    *pty,
                                   const char       *zname,
                                   cs_real_t         symtens[6])
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ANISO_SYM) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not anisotropic"
              " with a symmetric storage.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        6, /* dim */
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        symtens);

  pty->defs[new_id] = d;
  pty->get_eval_at_cell[new_id] = cs_xdef_eval_symtens_by_val;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_symtens_by_val;

  /* Set the state flag */

  pty->state_flag |= CS_FLAG_STATE_CELLWISE | CS_FLAG_STATE_STEADY;
  if (z_id == 0)
    pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value at the boundary for the given (anisotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *        Optimized case with a symmetric storage.
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]       tens     tensor to set given as an array of 6 values
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_aniso_sym_by_value(cs_property_t    *pty,
                                            const char       *zname,
                                            double            tens[6])
{
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if ((pty->type & CS_PROPERTY_ANISO_SYM) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property \"%s\" is not anisotropic"
              " with a symmetric storage.\n"
              " Please check your settings.", pty->name);

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE |
    CS_FLAG_STATE_STEADY;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          6,      /* dim */
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          tens); /* context */

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_property_t structure thanks to a time function for all
 *        cells associated to the zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           cells are considered)
 * \param[in]       func     pointer to a cs_time_func_t function
 * \param[in]       input    nullptr or pointer to a structure cast on-the-fly
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
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_time_func_context_t tfc
    = { .z_id = z_id, .func = func, .input = input, .free_input = nullptr };

  /* Default initialization */

  pty->get_eval_at_cell[new_id]    = nullptr;
  pty->get_eval_at_cell_cw[new_id] = cs_xdef_cw_eval_by_time_func;

  int  dim = 0;
  if (pty->type & CS_PROPERTY_ISO) {
    dim = 1;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_scalar_by_time_func;
  }
  else if (pty->type & CS_PROPERTY_ORTHO) {
    dim = 3;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_vector_by_time_func;
  }
  else if (pty->type & CS_PROPERTY_ANISO_SYM) {
    dim = 6;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_symtens_by_time_func;
  }
  else if (pty->type & CS_PROPERTY_ANISO) {
    dim = 9;
    pty->get_eval_at_cell[new_id] = cs_xdef_eval_tensor_by_time_func;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible type for the property \"%s\".",
              __func__, pty->name);

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_TIME_FUNCTION,
                                        dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &tfc);

  pty->defs[new_id] = d;

  /* Set the state flag for the property */

  pty->state_flag |= CS_FLAG_STATE_CELLWISE;
  if (z_id == 0)
    pty->state_flag |= CS_FLAG_STATE_UNIFORM;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure at the boundary thanks
 *        to a time function in a subdomain attached to a zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]       func     pointer to a cs_time_func_t function
 * \param[in]       input    nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_time_func(cs_property_t      *pty,
                                      const char         *zname,
                                      cs_time_func_t     *func,
                                      void               *input)
{
  if (func == nullptr)
    return nullptr;
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; /* metadata */
  cs_xdef_time_func_context_t tfc
    = { .z_id = z_id, .func = func, .input = input, .free_input = nullptr };

  int  dim = cs_property_get_dim(pty);

  if (dim == 0)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible type for the property \"%s\".",
              __func__, pty->name);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_TIME_FUNCTION,
                                          dim,
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          &tfc);

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_property_t structure thanks to an analytic function for
 *        all cells associated to the zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           cells are considered)
 * \param[in]       func     pointer to a cs_analytic_func_t function
 * \param[in]       input    nullptr or pointer to a structure cast on-the-fly
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
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  cs_flag_t  state_flag = 0, meta_flag = 0; /* metadata */
  int  z_id = cs_volume_zone_id_by_name(zname);
  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = func, .input = input, .free_input = nullptr };

  int  dim = cs_property_get_dim(pty);

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
 * \brief Define the value of a cs_property_t structure at the boundary thanks
 *        to a time function for all boundary faces associated to the zone
 *        named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]       func     pointer to a cs_analytic_func_t function
 * \param[in]       input    nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_analytic(cs_property_t        *pty,
                                     const char           *zname,
                                     cs_analytic_func_t   *func,
                                     void                 *input)
{
  if (func == nullptr)
    return nullptr;
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  int  dim = cs_property_get_dim(pty);
  cs_flag_t  state_flag = 0, meta_flag = 0; /* metadata */
  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = func, .input = input, .free_input = nullptr };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          dim,
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          &ac);

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure thanks to low-level
 *        functions specifying how to evaluate the property in each cell (with
 *        cell-wise structures or not). This definition applies to all cells
 *        associated to the zone named zname
 *
 * \param[in, out] pty           pointer to a cs_property_t structure
 * \param[in]      zname         name of the zone (if nullptr or "" then all
 * cells are selected) \param[in]      context       pointer to a structure (may
 * be nullptr) \param[in]      get_eval_at_cell      pointer to a function
 * \param[in]      get_eval_at_cell_cw   pointer to a function
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
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  def_id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0; /* metadata */

  int  dim = cs_property_get_dim(pty);

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
 * \brief Define a cs_property_t structure thanks to an array of values. If an
 *        advanced usage of the definition by array is needed, then call \ref
 *        cs_xdef_array_set_adjacency and/or \ref cs_xdef_array_set_sublist
 *
 * \param[in, out] pty           pointer to a cs_property_t structure
 * \param[in]      zname         name of the zone (if nullptr or "" then all
 * cells are selected) \param[in]      val_location  information to know where
 * are located values \param[in]      array         pointer to an array
 * \param[in]      is_owner      transfer the lifecycle to the cs_xdef_t struc.
 *                               (true or false)
 * \param[in]      full_length   if true, the size of "array" is allocated to
 *                               the total numbers of entities related to the
 *                               given location. If false, a new list is
 *                               allocated and filled with the related subset
 *                               indirection.
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_array(cs_property_t      *pty,
                         const char         *zname,
                         cs_flag_t           val_location,
                         cs_real_t          *array,
                         bool                is_owner,
                         bool                full_length)
{
  if (array == nullptr)
    return nullptr;
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  id = _add_new_def(pty);
  int  z_id = cs_volume_zone_id_by_name(zname);
  int  dim = cs_property_get_dim(pty);
  cs_flag_t  state_flag = 0; /* Will be updated during the creation */
  cs_flag_t  meta_flag = 0;  /* metadata */

  if (z_id == 0 && full_length == false) {
    full_length = true;
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  "%s: Invalid settings detected for property \"%s\"\n"
                  "    A full-length array is set since z_id=0.",
                  __func__, pty->name);
  }

  cs_xdef_array_context_t input = { .z_id           = z_id,
                                    .stride         = dim,
                                    .value_location = val_location,
                                    .is_owner       = is_owner,
                                    .full_length    = full_length,
                                    .values         = array,
                                    /* Optional parameters (not used) */
                                    .full2subset = nullptr,
                                    .n_list_elts = 0,
                                    .elt_ids     = nullptr,
                                    .adjacency   = nullptr };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &input);

  /* Build the indirection array if only a subset is used */

  if (!full_length)
    cs_xdef_array_build_full2subset(d);

  /* Set pointers */

  pty->defs[id] = d;

  if (dim == 1)
    pty->get_eval_at_cell[id] = cs_xdef_eval_scalar_at_cells_by_array;
  else
    pty->get_eval_at_cell[id] = cs_xdef_eval_nd_at_cells_by_array;
  pty->get_eval_at_cell_cw[id] = cs_xdef_cw_eval_by_array;

  if (cs_flag_test(val_location, cs_flag_primal_cell)     == false &&
      cs_flag_test(val_location, cs_flag_primal_vtx)      == false &&
      cs_flag_test(val_location, cs_flag_primal_edge_byc) == false &&
      cs_flag_test(val_location, cs_flag_dual_face_byc)   == false &&
      cs_flag_test(val_location, cs_flag_dual_cell_byc)   == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Property \"%s\". Case not available.\n",
              __func__, pty->name);

  /* Set the state flag */

  if (cs_flag_test(val_location, cs_flag_primal_cell))
    pty->state_flag |= CS_FLAG_STATE_CELLWISE;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the values of a property at the boundary thanks to an array.
 *        If an advanced usage of the definition by array is needed, then call
 *        \ref cs_xdef_array_set_adjacency and/or \ref
 *        cs_xdef_array_set_sublist
 *
 * \param[in, out] pty          pointer to a cs_property_t structure
 * \param[in]      zname        nullptr or name of the boundary zone
 * \param[in]      val_loc      information to know where are located values
 * \param[in]      array        pointer to an array
 * \param[in]      is_owner     transfer the lifecycle to the cs_xdef_t struct.
 *                              (true or false)
 * \param[in]      full_length  if true, the size of "array" should be allocated
 *                              to the total numbers of entities related to the
 *                              given location. If false, a new list is
 *                              allocated and filled with the related subset
 *                              indirection.
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_array(cs_property_t      *pty,
                                  const char         *zname,
                                  cs_flag_t           val_loc,
                                  cs_real_t          *array,
                                  bool                is_owner,
                                  bool                full_length)
{
  if (array == nullptr)
    return nullptr;
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_b_def(pty);
  int  z_id = cs_boundary_zone_id_by_name(zname);
  int  dim = cs_property_get_dim(pty);
  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0;

  if (z_id == 0 && full_length == false) {
    full_length = true;
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  "%s: Invalid settings detected for property \"%s\"\n"
                  "    A full-length array is set since z_id=0.",
                  __func__, pty->name);
  }

  cs_xdef_array_context_t input = { .z_id           = z_id,
                                    .stride         = dim,
                                    .value_location = val_loc,
                                    .is_owner       = is_owner,
                                    .full_length    = full_length,
                                    .values         = array,
                                    /* Optional parameters (not used) */
                                    .full2subset = nullptr,
                                    .n_list_elts = 0,
                                    .elt_ids     = nullptr,
                                    .adjacency   = nullptr };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          dim,
                                          z_id, /* zone_id */
                                          state_flag,
                                          meta_flag,
                                          &input);

  /* Build the indirection array if only a subset is used */

  if (!full_length)
    cs_xdef_array_build_full2subset(d);

  pty->b_defs[new_id] = d;

  if (cs_flag_test(val_loc, cs_flag_primal_face)   == false &&
      cs_flag_test(val_loc, cs_flag_primal_vtx)    == false &&
      cs_flag_test(val_loc, cs_flag_boundary_face) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Property \"%s\". Case not available.\n",
              __func__, pty->name);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a field structure. One
 *         assumes that all cells are defined using this array.
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_field(cs_property_t    *pty,
                         cs_field_t       *field)
{
  if (field == nullptr)
    return nullptr;
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  id = _add_new_def(pty);
  int  dim = cs_property_get_dim(pty);
  cs_flag_t  state_flag = CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; /* metadata */

  assert(dim == field->dim);
  assert(id == 0); /* z_id = 0 since all the support is selected in this case */

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

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_FIELD,
                                        dim,
                                        0, /* zone_id */
                                        state_flag,
                                        meta_flag,
                                        field);

  pty->defs[id] = d;
  pty->get_eval_at_cell[id] = cs_xdef_eval_cell_by_field;
  pty->get_eval_at_cell_cw[id] = cs_xdef_cw_eval_by_field;

  /* Set the state flag */

  pty->state_flag |= CS_FLAG_STATE_CELLWISE;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the values of a property at the boundary thanks to a field
 *        structure. One assumes that all boundary faces are defined using
 *        this array.
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_field(cs_property_t    *pty,
                                  cs_field_t       *field)
{
  if (field == nullptr)
    return nullptr;
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = _add_new_b_def(pty);
  int  dim = cs_property_get_dim(pty);

  /* Only one definition if one uses a definition by field. No zone name is
     given since one assumes that z_id = 0 */

  if (pty->n_definitions > 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: When a definition by field is requested, the max. number"
              " of zones to consider should be equal to 1.\n"
              " Current value is %d for property \"%s\".\n"
              " Please check your settings.",
              __func__, pty->n_definitions, pty->name);

  const cs_zone_t  *z = cs_boundary_zone_by_id(0);

  assert(dim == field->dim);
  if (field->location_id != z->location_id)
    bft_error(__FILE__, __LINE__, 0,
              " Property defined by field requests that the field location"
              " is supported by cells\n"
              " Property %s\n", pty->name);

  cs_flag_t  state_flag = CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = 0; /* metadata */

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_FIELD,
                                          dim,
                                          0, /* zone_id */
                                          state_flag,
                                          meta_flag,
                                          field);

  pty->b_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the values of the property at cells from the given
 *        definition. According to the parameter "dense_ouput", the "eval"
 *        array should be allocated with a size equal to pty->dim*n_cells
 *        (where "dim" depends on the type of property to handle) when no dense
 *        ouput is requested. Otherwise, an allocation size equal to pty->dim *
 *        the number of cells associated to the definition "def" is enough.
 *
 *        No scaling is applied to the value. This should be done with a higher
 *        level function like \ref cs_property_eval_at_cells or
 *        \ref cs_property_get_cell_tensor
 *
 * \param[in]      pty           pointer to a property structure
 * \param[in]      def_id        id associated to the definition
 * \param[in]      dense_output  true/false
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in, out] eval          array storing the result of the evaluations
 */
/*----------------------------------------------------------------------------*/

void
cs_property_evaluate_def(const cs_property_t    *pty,
                         int                     def_id,
                         bool                    dense_output,
                         double                  t_eval,
                         cs_real_t              *eval)
{
  if (pty == nullptr)
    return;
  assert(def_id > -1 && def_id < pty->n_definitions);

  const cs_xdef_t  *def = pty->defs[def_id];
  if (def == nullptr)
    return;

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

  if (def->type == CS_XDEF_BY_SUB_DEFINITIONS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of definition. Property \"%s\"; Zone \"%s\".\n",
              __func__, pty->name, z->name);

  if (def->z_id != 0) /* Not the full support */
    pty->get_eval_at_cell[def_id](z->n_elts,
                                  z->elt_ids,
                                  dense_output,
                                  cs_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  eval);

  else /* All elements are selected: elt_ids = nullptr */
    pty->get_eval_at_cell[def_id](z->n_elts,
                                  nullptr,
                                  dense_output,
                                  cs_mesh,
                                  cs_cdo_connect,
                                  cs_cdo_quant,
                                  t_eval,
                                  def->context,
                                  eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the values of a property at boundary faces from the given
 *        boundary definition. If a dense ouput is not requested, then the size
 *        of the resulting array should be allocated at least to pty->dim *
 *        n_b_faces. Otherwise, n_b_faces can be replaced by the number of
 *        boundary faces associated to the current definition.
 *
 *        No scaling is applied to the value. This should be done with a higher
 *        level function like \ref cs_property_eval_at_boundary_faces
 *
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      def_id        id associated to the definition
 * \param[in]      dense_output  true/false
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in, out] array         array storing the result of the evaluation(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_evaluate_boundary_def(const cs_property_t  *pty,
                                  int                   def_id,
                                  bool                  dense_output,
                                  double                t_eval,
                                  cs_real_t            *array)
{
  if (pty == nullptr)
    return;
  assert(def_id > -1 && def_id < pty->n_b_definitions);

  const cs_xdef_t  *def = pty->b_defs[def_id];
  if (def == nullptr)
    return;

  assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);
  int  pty_dim = cs_property_get_dim(pty);

  const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);

  cs_lnum_t  n_elts = z->n_elts;
  const cs_lnum_t  *elt_ids;

  if (def->z_id != 0) /* Not the full support */
    elt_ids = z->elt_ids;
  else
    elt_ids = nullptr;

  switch (def->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    /* ---------------------------- */
    {
      cs_xdef_analytic_context_t *cx
        = static_cast<cs_xdef_analytic_context_t *>(def->context);
      assert(cx != nullptr);

      switch (def->qtype) {

      case CS_QUADRATURE_BARY_SUBDIV:
      case CS_QUADRATURE_HIGHER:
      case CS_QUADRATURE_HIGHEST:
        cs_xdef_eval_avg_at_b_faces_by_analytic(n_elts, elt_ids,
                                                dense_output,
                                                cs_mesh,
                                                cs_cdo_connect,
                                                cs_cdo_quant,
                                                t_eval,
                                                cx,
                                                def->qtype,
                                                pty_dim,
                                                array);
        break;

      default:
        cs_xdef_eval_at_b_faces_by_analytic(n_elts, elt_ids,
                                            dense_output,
                                            cs_mesh,
                                            cs_cdo_connect,
                                            cs_cdo_quant,
                                            t_eval,
                                            def->context,
                                            array);
        break;
      }
    }
    break;

  case CS_XDEF_BY_ARRAY:
    /* ---------------- */
    {
      cs_xdef_array_context_t *cx
        = static_cast<cs_xdef_array_context_t *>(def->context);
      assert(cx->stride == pty_dim);

      if (cs_flag_test(cx->value_location, cs_flag_primal_face)) {

        if (elt_ids == nullptr) /* All the boundary is considered */
          cs_array_real_copy(pty_dim*cs_mesh->n_b_faces, cx->values, array);

        else { /* Only a part of the boundary is considered */

          if (cx->full_length) {

            if (dense_output)
              cs_array_real_copy_subset(n_elts, pty_dim, elt_ids,
                                        CS_ARRAY_SUBSET_IN,
                                        cx->values,
                                        array);
            else
              cs_array_real_copy_subset(n_elts, pty_dim, elt_ids,
                                        CS_ARRAY_SUBSET_INOUT,
                                        cx->values,
                                        array);

          }
          else {

            if (dense_output)
              cs_array_real_copy(pty_dim*n_elts, cx->values, array);
            else
              cs_array_real_copy_subset(n_elts, pty_dim, elt_ids,
                                        CS_ARRAY_SUBSET_OUT,
                                        cx->values,
                                        array);

          }

        } /* Only a part of the boundary is defined */

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid location. Property \"%s\"; Zone \"%s\".\n",
                  __func__, pty->name, z->name);
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    /* ----------------------- */
    {
      const cs_xdef_dof_context_t *cx
        = static_cast<cs_xdef_dof_context_t *>(def->context);
      assert(cx != nullptr);

      if (cs_flag_test(cx->dof_location, cs_flag_primal_face))
        cx->func(n_elts, elt_ids, dense_output, cx->input, array);
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid location. Property \"%s\"; Zone \"%s\".\n",
                  __func__, pty->name, z->name);
    }
    break;

  case CS_XDEF_BY_FIELD:
    /* ---------------- */
    {
      cs_field_t  *field = (cs_field_t *)def->context;

      if (cs_mesh_location_get_type(field->location_id) !=
          CS_MESH_LOCATION_BOUNDARY_FACES)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid field location. Property \"%s\"; Zone \"%s\".\n",
                  __func__, pty->name, z->name);

      cs_array_real_copy(field->dim*cs_mesh->n_b_faces, field->val, array);
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    /* ------------------------ */
    {
      const cs_xdef_time_func_context_t *tfc
        = static_cast<cs_xdef_time_func_context_t *>(def->context);
      assert(tfc != nullptr);

      cs_real_t  ref_val[9];
      tfc->func(t_eval, tfc->input, ref_val);

      if (dense_output)
        _assign_ref_value(pty_dim, n_elts, nullptr, ref_val, array);
      else
        _assign_ref_value(pty_dim, n_elts, elt_ids, ref_val, array);
    }
    break;

  case CS_XDEF_BY_VALUE:
    /* ---------------- */
    {
      const cs_real_t *ref_val = static_cast<const cs_real_t *>(def->context);

      if (dense_output)
        _assign_ref_value(pty_dim, n_elts, nullptr, ref_val, array);
      else
        _assign_ref_value(pty_dim, n_elts, elt_ids, ref_val, array);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid definition type. Property \"%s\"; Zone \"%s\".\n",
              __func__, pty->name, z->name);

  } /* Type of definition */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of the property at each cell. Store the evaluation
 *        in the given array.
 *
 * \param[in]      t_eval       physical time at which one evaluates the term
 * \param[in]      pty          pointer to a cs_property_t structure
 * \param[out]     pty_stride   = 0 if uniform, =1 otherwise
 * \param[in, out] p_pty_vals   pointer to an array of values. Allocated if not
 *                              The size of the allocation depends on the value
 *                              of the pty_stride
 */
/*----------------------------------------------------------------------------*/

void
cs_property_iso_get_cell_values(cs_real_t               t_eval,
                                const cs_property_t    *pty,
                                int                    *pty_stride,
                                cs_real_t             **p_pty_vals)
{
  if (pty == nullptr)
    return;
  assert(pty->type & CS_PROPERTY_ISO);

  bool        allocate = (*p_pty_vals == nullptr) ? true : false;
  cs_real_t  *values = *p_pty_vals;

  if (cs_property_is_uniform(pty)) {

    *pty_stride = 0;
    if (allocate)
      BFT_MALLOC(values, 1, cs_real_t);

    /* Evaluation at c_id = 0. One assumes that there is at least one cell per
       MPI rank */

    values[0] = cs_property_get_cell_value(0, t_eval, pty);

    /* scaling is performed if requested inside the previous call */

  }
  else {

    *pty_stride = 1;
    if (allocate)
      BFT_MALLOC(values, cs_cdo_quant->n_cells, cs_real_t);

    cs_property_eval_at_cells(t_eval, pty, values);

    /* scaling is performed if requested inside the previous call */

  }

  /* Return the pointer to values */

  *p_pty_vals = values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of the property at each cell. Store the
 *        evaluation in the given array.
 *
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      pty      pointer to a cs_property_t structure
 * \param[in, out] array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_cells(cs_real_t               t_eval,
                          const cs_property_t    *pty,
                          cs_real_t              *array)
{
  if (pty == nullptr)
    return;
  assert(array != nullptr);

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_cells = quant->n_cells;

  double  scaling_factor =
    (pty->type & CS_PROPERTY_SCALED) ? pty->scaling_factor : 1.0;

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    assert(pty->related_properties != nullptr);
    const cs_property_t  *a = pty->related_properties[0];
    const cs_property_t  *b = pty->related_properties[1];

    if (a->type & CS_PROPERTY_SCALED)
      scaling_factor *= a->scaling_factor;
    if (b->type & CS_PROPERTY_SCALED)
      scaling_factor *= b->scaling_factor;

    cs_real_t *tmp_val = nullptr;
    BFT_MALLOC(tmp_val, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, tmp_val);

    if (pty->type & CS_PROPERTY_ISO) {

      /* 1. Evaluate the property A */

      for (int def_id = 0; def_id < a->n_definitions; def_id++)
        cs_property_evaluate_def(a,
                                 def_id,
                                 false, /* dense output */
                                 t_eval,
                                 tmp_val);

      /* 2. Evaluate the property B */

      for (int def_id = 0; def_id < b->n_definitions; def_id++)
        cs_property_evaluate_def(b,
                                 def_id,
                                 false, /* dense output */
                                 t_eval,
                                 array);

      /* 3. Operate the product of the two evaluations */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        array[i] *= tmp_val[i];

    }
    else {

      if (a->type & CS_PROPERTY_ISO) {

        /* 1. Evaluate the property A */

        for (int def_id = 0; def_id < a->n_definitions; def_id++)
          cs_property_evaluate_def(a,
                                   def_id,
                                   false, /* dense output */
                                   t_eval,
                                   tmp_val);

        /* 2. Evaluate the property B */

        int  b_dim = cs_property_get_dim(b);

        for (int def_id = 0; def_id < b->n_definitions; def_id++)
          cs_property_evaluate_def(b,
                                   def_id,
                                   false, /* dense output */
                                   t_eval,
                                   array);

        /* 3. Operate the product of the two evaluations */

        for (cs_lnum_t i = 0; i < n_cells; i++) {
          const cs_real_t  acoef = tmp_val[i];
          cs_real_t  *_a = array + b_dim*i;
          for (int k = 0; k < b_dim; k++)
            _a[k] *= acoef;
        }

      }
      else if (b->type & CS_PROPERTY_ISO) {

        /* 1. Evaluate the property B */

        for (int def_id = 0; def_id < b->n_definitions; def_id++)
          cs_property_evaluate_def(b,
                                   def_id,
                                   false, /* dense output */
                                   t_eval,
                                   tmp_val);

        /* 2. Evaluate the property A */

        int  a_dim = cs_property_get_dim(a);

        for (int def_id = 0; def_id < a->n_definitions; def_id++)
          cs_property_evaluate_def(a,
                                   def_id,
                                   false, /* dense output */
                                   t_eval,
                                   array);

        /* 3. Operate the product of the two evaluations */

        for (cs_lnum_t i = 0; i < n_cells; i++) {
          const cs_real_t  bcoef = tmp_val[i];
          cs_real_t  *_a = array + a_dim*i;
          for (int k = 0; k < a_dim; k++)
            _a[k] *= bcoef;
        }

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Property \"%s\". Case not handled yet.\n",
                  __func__, pty->name);

    } /* Either a or b is an isotropic property */

    BFT_FREE(tmp_val);

  }
  else { /* Simple case: One has to evaluate the property */

    if ((pty->type & CS_PROPERTY_ISO) && cs_property_is_constant(pty))
      cs_array_real_set_scalar(n_cells, pty->ref_value, array);

    else {

      for (int def_id = 0; def_id < pty->n_definitions; def_id++)
        cs_property_evaluate_def(pty,
                                 def_id,
                                 false, /* dense output */
                                 t_eval,
                                 array);

    }

  } /* Not defined as the product of two existing properties */

  /* Apply a scaling factor is requested */

  if (fabs(scaling_factor - 1.0) > 10*FLT_MIN)
    cs_array_real_scale(
      n_cells, cs_property_get_dim(pty), nullptr, scaling_factor, array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of the property at each boundary face. Store the
 *        result of the evaluation in the given array.
 *
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      pty      pointer to a cs_property_t structure
 * \param[in, out] array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_boundary_faces(cs_real_t               t_eval,
                                   const cs_property_t    *pty,
                                   cs_real_t              *array)
{
  if (pty == nullptr)
    return;
  assert(array != nullptr);

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_b_faces = quant->n_b_faces;

  double  scaling_factor =
    (pty->type & CS_PROPERTY_SCALED) ? pty->scaling_factor : 1.0;

  if (pty->n_b_definitions == 0) {

    /* One uses the definition in the cell associated to each boundary face to
       get the result */

    if (pty->type & CS_PROPERTY_BY_PRODUCT) { /* P = A * B */

      assert(pty->related_properties != nullptr);
      const cs_property_t  *a = pty->related_properties[0];
      const cs_property_t  *b = pty->related_properties[1];

      if (a->type & CS_PROPERTY_SCALED)
        scaling_factor *= a->scaling_factor;
      if (b->type & CS_PROPERTY_SCALED)
        scaling_factor *= b->scaling_factor;

      cs_lnum_t *a_def_idx = nullptr, *a_cell_ids = nullptr,
                *a_bf_ids  = nullptr;
      cs_lnum_t *b_def_idx = nullptr, *b_cell_ids = nullptr,
                *b_bf_ids = nullptr;

      _build_def_idx_from_bdy_selection(a, &a_def_idx, &a_cell_ids, &a_bf_ids);
      _build_def_idx_from_bdy_selection(b, &b_def_idx, &b_cell_ids, &b_bf_ids);

      int  dim_a = cs_property_get_dim(a);
      int  dim_b = cs_property_get_dim(b);
      int  dim_ab = CS_MAX(dim_a, dim_b);

      cs_real_t *tmp_val = nullptr;
      BFT_MALLOC(tmp_val, n_b_faces*dim_ab, cs_real_t);
      cs_array_real_fill_zero(n_b_faces*dim_ab, tmp_val);

      if (pty->type & CS_PROPERTY_ISO) { /* A and B are isotropic */

        /* 1. Evaluates the property A */

        if (a->n_definitions == 1)
          _evaluate_property_at_boundary_from_cells(a, a_def_idx, a_cell_ids,
                                                    t_eval,
                                                    array);

        else { /* Several definitions for the property A */

          _evaluate_property_at_boundary_from_cells(a, a_def_idx, a_cell_ids,
                                                    t_eval,
                                                    tmp_val);

          /* Apply the indirection to fill the array to return */

          cs_array_real_copy_subset(n_b_faces, 1, a_bf_ids,
                                    CS_ARRAY_SUBSET_IN,
                                    tmp_val,
                                    array);

        }

        /* 2. Evaluates the property B and operates the product */

        _evaluate_property_at_boundary_from_cells(b, b_def_idx, b_cell_ids,
                                                  t_eval,
                                                  tmp_val);

        if (b->n_definitions == 1) {

          for (cs_lnum_t i = 0; i < n_b_faces; i++)
            array[i] *= tmp_val[i];

        }
        else { /* Several definitions for the property B */

          for (cs_lnum_t i = 0; i < n_b_faces; i++)
            array[i] *= tmp_val[b_bf_ids[i]];

        }

      }
      else {

        if (a->type & CS_PROPERTY_ISO) {

          /* 1. Evaluates the property B (not isotropic) */

          if (b->n_definitions == 1)
            _evaluate_property_at_boundary_from_cells(b, b_def_idx, b_cell_ids,
                                                      t_eval,
                                                      array);

          else { /* Several definitions for the property B */

            _evaluate_property_at_boundary_from_cells(b, b_def_idx, b_cell_ids,
                                                      t_eval,
                                                      tmp_val);

            /* Apply the indirection to fill the array to return */

            cs_array_real_copy_subset(n_b_faces, dim_b, b_bf_ids,
                                      CS_ARRAY_SUBSET_IN,
                                      tmp_val,
                                      array);

          }

          /* 2. Evaluates the property A (isotropic) and operates the product */

          _evaluate_property_at_boundary_from_cells(a, a_def_idx, a_cell_ids,
                                                    t_eval,
                                                    tmp_val);

          if (a->n_definitions == 1) {

            for (cs_lnum_t i = 0; i < n_b_faces; i++) {

              const cs_real_t coef = tmp_val[i]; /* value of the A property */
              cs_real_t  *_array = array + i*dim_b;
              for (int k = 0; k < dim_b; k++)
                _array[k] *= coef;

            }

          }
          else { /* Several definitions for the property A */

            for (cs_lnum_t i = 0; i < n_b_faces; i++) {

              const cs_real_t coef = tmp_val[a_bf_ids[i]];
              cs_real_t  *_array = array + i*dim_b;
              for (int k = 0; k < dim_b; k++)
                _array[k] *= coef;

            }

          }

        }
        else if (b->type & CS_PROPERTY_ISO) {

          /* 1. Evaluates the property A (not isotropic) */

          if (a->n_definitions == 1)
            _evaluate_property_at_boundary_from_cells(a, a_def_idx, a_cell_ids,
                                                      t_eval,
                                                      array);

          else { /* Several definitions for the property A */

            _evaluate_property_at_boundary_from_cells(a, a_def_idx, a_cell_ids,
                                                      t_eval,
                                                      tmp_val);

            /* Apply the indirection to fill the array to return */

            cs_array_real_copy_subset(n_b_faces, dim_a, a_bf_ids,
                                      CS_ARRAY_SUBSET_IN,
                                      tmp_val,
                                      array);

          }

          /* 2. Evaluates the property B (isotropic) and operates the product */

          _evaluate_property_at_boundary_from_cells(b, b_def_idx, b_cell_ids,
                                                    t_eval,
                                                    tmp_val);

          if (b->n_definitions == 1) {

            for (cs_lnum_t i = 0; i < n_b_faces; i++) {

              const cs_real_t coef = tmp_val[i]; /* value of the B property */
              cs_real_t  *_array = array + i*dim_a;
              for (int k = 0; k < dim_a; k++)
                _array[k] *= coef;

            }

          }
          else { /* Several definitions for the property B */

            for (cs_lnum_t i = 0; i < n_b_faces; i++) {

              const cs_real_t coef = tmp_val[b_bf_ids[i]];
              cs_real_t  *_array = array + i*dim_a;
              for (int k = 0; k < dim_a; k++)
                _array[k] *= coef;

            }

          }

        }
        else
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Property \"%s\". Case not handled yet.\n",
                    __func__, pty->name);

      } /* Either a or b is an isotropic property */

      BFT_FREE(tmp_val);

    }
    else { /* Simple case: One has to evaluate the property */

      if ((pty->type & CS_PROPERTY_ISO) && cs_property_is_constant(pty))
        cs_array_real_set_scalar(n_b_faces, pty->ref_value, array);

      else { /* One property but less simple case */

        cs_lnum_t *def_idx = nullptr, *cell_ids = nullptr, *bf_ids = nullptr;

        if (pty->n_definitions == 1)
          _evaluate_property_at_boundary_from_cells(pty, def_idx, cell_ids,
                                                    t_eval,
                                                    array);

        else {

          _build_def_idx_from_bdy_selection(pty, &def_idx, &cell_ids, &bf_ids);

          int  dim = cs_property_get_dim(pty);
          cs_real_t *tmp_val = nullptr;
          BFT_MALLOC(tmp_val, n_b_faces*dim, cs_real_t);
          cs_array_real_fill_zero(n_b_faces*dim, tmp_val);

          _evaluate_property_at_boundary_from_cells(pty, def_idx, cell_ids,
                                                    t_eval,
                                                    tmp_val);

          /* Apply the indirection to fill the array to return */

          cs_array_real_copy_subset(n_b_faces, dim, bf_ids,
                                    CS_ARRAY_SUBSET_IN,
                                    tmp_val,
                                    array);

        }

      } /* Not isotropic and constant at the same time */

    } /* Not defined as the product of two existing properties */

  }
  else { /* There is at least one boundary definition */

    assert(pty->b_defs != nullptr);
    for (int def_id = 0; def_id < pty->n_b_definitions; def_id++)
      cs_property_evaluate_boundary_def(pty,
                                        def_id,
                                        false, /* dense ouput ? */
                                        t_eval,
                                        array);

  } /* n_b_definitions > 0 */

  /* Apply a scaling factor is requested */

  if (fabs(scaling_factor - 1.0) > 10*FLT_MIN)
    cs_array_real_scale(
      n_b_faces, cs_property_get_dim(pty), nullptr, scaling_factor, array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of the tensor attached a property at the cell
 *        center
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
  if (pty == nullptr)
    return;

  /* Initialize extra-diag. values of the tensor */

  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    _get_cell_tensor_by_property_product(c_id, t_eval, pty, tensor);

    const cs_property_t  *pty_a = pty->related_properties[0];
    if (pty_a->type & CS_PROPERTY_SCALED) {
      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          tensor[ki][kj] *= pty_a->scaling_factor;
    }

    const cs_property_t  *pty_b = pty->related_properties[1];
    if (pty_b->type & CS_PROPERTY_SCALED) {
      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          tensor[ki][kj] *= pty_b->scaling_factor;
    }

  }
  else
    _get_cell_tensor(c_id, t_eval, pty, tensor);

  if (pty->type & CS_PROPERTY_SCALED) {
    for (int ki = 0; ki < 3; ki++)
      for (int kj = 0; kj < 3; kj++)
        tensor[ki][kj] *= pty->scaling_factor;
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

  if (pty == nullptr)
    return result;

  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", __func__, pty->name);

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    assert(pty->related_properties != nullptr);

    const cs_property_t  *pty_a = pty->related_properties[0];
    result = _get_cell_value(c_id, t_eval, pty_a);
    if (pty_a->type & CS_PROPERTY_SCALED)
      result *= pty_a->scaling_factor;

    const cs_property_t  *pty_b = pty->related_properties[1];
    result *= _get_cell_value(c_id, t_eval, pty_b);
    if (pty_b->type & CS_PROPERTY_SCALED)
      result *= pty_b->scaling_factor;

  }
  else {

    if (cs_property_is_constant(pty))
      result = pty->ref_value;
    else
      result = _get_cell_value(c_id, t_eval, pty);

  }

  if (pty->type & CS_PROPERTY_SCALED)
    result *= pty->scaling_factor;

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of the tensor attached to a property at the cell
 *        center
 *        Version using a cs_cell_mesh_t structure
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
  if (pty == nullptr)
    return;

  /* Initialize extra-diag. values of the tensor */

  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    _tensor_in_cell_by_property_product(cm, pty, t_eval, tensor);

    const cs_property_t  *pty_a = pty->related_properties[0];
    if (pty_a->type & CS_PROPERTY_SCALED) {
      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          tensor[ki][kj] *= pty_a->scaling_factor;
    }

    const cs_property_t  *pty_b = pty->related_properties[1];
    if (pty_b->type & CS_PROPERTY_SCALED) {
      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          tensor[ki][kj] *= pty_b->scaling_factor;
    }

  }
  else
    _tensor_in_cell(cm, pty, t_eval, tensor);

  if (pty->type & CS_PROPERTY_SCALED) {
    for (int ki = 0; ki < 3; ki++)
      for (int kj = 0; kj < 3; kj++)
        tensor[ki][kj] *= pty->scaling_factor;
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
 * \brief Compute the value of a property at the cell center
 *        Version using a cs_cell_mesh_t structure
 *
 * \param[in] cm        pointer to a cs_cell_mesh_t structure
 * \param[in] pty       pointer to a cs_property_t structure
 * \param[in] t_eval    physical time at which one evaluates the term
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

  if (pty == nullptr)
    return result;

  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    assert(pty->related_properties != nullptr);

    const cs_property_t  *pty_a = pty->related_properties[0];
    result = _value_in_cell(cm, pty_a, t_eval);
    if (pty_a->type & CS_PROPERTY_SCALED)
      result *= pty_a->scaling_factor;

    const cs_property_t  *pty_b = pty->related_properties[1];
    result *= _value_in_cell(cm, pty_b, t_eval);
    if (pty_b->type & CS_PROPERTY_SCALED)
      result *= pty_b->scaling_factor;

  }
  else {

    if (cs_property_is_constant(pty))
      result = pty->ref_value;
    else
      result = _value_in_cell(cm, pty, t_eval);

  }

  if (pty->type & CS_PROPERTY_SCALED)
    result *= pty->scaling_factor;

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of an isotropic property in each portion of dual
 *        cell in a (primal) cell. This relies on the c2v connectivity.
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty       pointer to a cs_property_t structure
 * \param[in]      t_eval    physical time at which one evaluates the term
 * \param[in, out] eval      array of values storing the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_property_c2v_values(const cs_cell_mesh_t   *cm,
                       const cs_property_t    *pty,
                       cs_real_t               t_eval,
                       cs_real_t              *eval)
{
  if (pty == nullptr)
    return;

  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  assert(cs_property_is_subcell(pty));
  assert(eval != nullptr);

  const cs_xdef_t *def = pty->defs[_get_cell_def_id(cm->c_id, pty)];

  switch(def->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_context_t  *cx =
        (cs_xdef_analytic_context_t *)def->context;
      assert(cx != nullptr);

      /* Evaluate the function for this time inside each portion of dual cell */

      for (int i = 0; i < cm->n_vc; i++) {

        double xvc[3] = { 0.5*cm->xc[0] + 0.5*cm->xv[3*i],
                          0.5*cm->xc[1] + 0.5*cm->xv[3*i+1],
                          0.5*cm->xc[2] + 0.5*cm->xv[3*i+2] };

        cx->func(t_eval, 1, nullptr, xvc, true, cx->input, eval + i);
      }
    }
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)def->context;

      assert(cx->stride == 1);

      if (cx->value_location == cs_flag_primal_vtx ||
          cx->value_location == cs_flag_dual_cell) {

        for (int i = 0; i < cm->n_vc; i++)
          eval[i] = cx->values[cm->v_ids[i]];

      }
      else if (cx->value_location == cs_flag_primal_cell) {

        const cs_real_t val_c = cx->values[cm->c_id];
        for (int i = 0; i < cm->n_vc; i++)
          eval[i] = val_c;

      }
      else if (cx->value_location == cs_flag_dual_cell_byc) {

        const cs_adjacency_t  *adj = cx->adjacency;
        cs_real_t  *_val = cx->values + adj->idx[cm->c_id];
        memcpy(eval, _val, cm->n_vc*sizeof(cs_real_t));

      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid location for array.",
                  __func__);
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *fld = (cs_field_t *)def->context;

      assert(fld != nullptr);
      assert(fld->dim == 1);

      if (fld->location_id == cs_mesh_location_get_id_by_name(N_("cells"))) {

        const cs_real_t val_c = fld->val[cm->c_id];
        for (int i = 0; i < cm->n_vc; i++)
          eval[i] = val_c;

      }
      else if (fld->location_id ==
               cs_mesh_location_get_id_by_name(N_("vertices"))) {

        for (int i = 0; i < cm->n_vc; i++)
          eval[i] = fld->val[cm->v_ids[i]];

      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid location for field.",
                  __func__);
    }
    break;

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->context;

      for (int i = 0; i < cm->n_vc; i++)
        eval[i] = constant_val[0];
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid definition.", __func__);

  } /* Type of definition */

  if (pty->type & CS_PROPERTY_SCALED)
    for (int i = 0; i < cm->n_vc; i++)
      eval[i] *= pty->scaling_factor;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of a property at the face center
 *        by analytic if possible or by averaging value at cells center
 *
 * \param[in] f_id     id of the current face
 * \param[in] t_eval   physical time at which one evaluates the term
 * \param[in] pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given face
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_get_face_value(cs_lnum_t            f_id,
                           cs_real_t            t_eval,
                           const cs_property_t *pty)
{
  cs_real_t result = 0;

  if (pty == nullptr)
    return result;

  if ((pty->type & CS_PROPERTY_ISO) == 0)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid type of property for this function.\n"
              " Property %s has to be isotropic.",
              __func__,
              pty->name);

  if (pty->type & CS_PROPERTY_BY_PRODUCT) {

    assert(pty->related_properties != nullptr);

    const cs_property_t *pty_a = pty->related_properties[0];
    result                     = _get_face_value(f_id, t_eval, pty_a);
    if (pty_a->type & CS_PROPERTY_SCALED)
      result *= pty_a->scaling_factor;

    const cs_property_t *pty_b = pty->related_properties[1];
    result *= _get_face_value(f_id, t_eval, pty_b);
    if (pty_b->type & CS_PROPERTY_SCALED)
      result *= pty_b->scaling_factor;
  }
  else {

    if (cs_property_is_constant(pty))
      result = pty->ref_value;
    else
      result = _get_face_value(f_id, t_eval, pty);
  }

  if (pty->type & CS_PROPERTY_SCALED)
    result *= pty->scaling_factor;

  return result;
};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Fourier number in each cell
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
  assert(fourier != nullptr); /* Sanity check */
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
 * \brief Print a summary of the settings for all defined cs_property_t
 *        structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_log_setup(void)
{
  if (_n_properties == 0)
    return;
  assert(_properties != nullptr);

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the definition of properties\n");
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h1);

  char  prefix[256];

  for (int i = 0; i < _n_properties; i++) {

    bool  is_uniform = false, is_steady = false;
    const cs_property_t  *pty = _properties[i];

    if (pty == nullptr)
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
    else if (pty->type & CS_PROPERTY_ANISO_SYM)
      cs_log_printf(CS_LOG_SETUP,
                    "  * %s | Type: anisotropic with a symmetric storage",
                    pty->name);
    else if (pty->type & CS_PROPERTY_ANISO)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Type: anisotropic", pty->name);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: Invalid type of property."),
                __func__);

    if (pty->type & CS_PROPERTY_BY_PRODUCT)
      cs_log_printf(CS_LOG_SETUP, " | by product\n");
    else
      cs_log_printf(CS_LOG_SETUP, "\n");

    cs_log_printf(CS_LOG_SETUP, "  * %s | Subcell definition %s\n",
                  pty->name, cs_base_strtf(cs_property_is_subcell(pty)));

    cs_log_printf(CS_LOG_SETUP, "  * %s | Number of definitions: %d\n\n",
                  pty->name, pty->n_definitions);

    for (int j = 0; j < pty->n_definitions; j++) {
      sprintf(prefix, "        Definition %3d", j);
      cs_xdef_log_setup(prefix, pty->defs[j]);
    }

  } /* Loop on properties */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
