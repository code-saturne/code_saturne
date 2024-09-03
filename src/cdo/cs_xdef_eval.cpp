/*============================================================================
 * Manage the (generic) evaluation of extended definitions
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_defs.h"
#include "cs_field.h"
#include "cs_mesh_location.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_xdef_eval.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */

#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a scalar-valued quantity for a list of elements
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_scalar_by_val(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         dense_output,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *context,
                           cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  const cs_real_t  *constant_val = (cs_real_t *)context;
  assert(eval != nullptr && constant_val != nullptr);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_scalar_on_subset(n_elts, elt_ids, constant_val[0], eval);
  else
    cs_array_real_set_scalar(n_elts, constant_val[0], eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a vector-valued quantity for a list of elements
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_vector_by_val(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         dense_output,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *context,
                           cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  const cs_real_t  *constant_val = (cs_real_t *)context;
  assert(eval != nullptr && constant_val != nullptr);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_vector_on_subset(n_elts, elt_ids, constant_val, eval);
  else
    cs_array_real_set_vector(n_elts, constant_val, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a tensor-valued quantity for a list of elements with
 *         symmetric storage.
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_symtens_by_val(cs_lnum_t                    n_elts,
                            const cs_lnum_t             *elt_ids,
                            bool                         dense_output,
                            const cs_mesh_t             *mesh,
                            const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *quant,
                            cs_real_t                    time_eval,
                            void                        *context,
                            cs_real_t                   *eval)
{
  CS_UNUSED(quant);
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  const cs_real_t  *constant_val = (const cs_real_t *)context;
  assert(eval != nullptr && constant_val != nullptr);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_value_on_subset(n_elts, 6, elt_ids, constant_val, eval);
  else
    cs_array_real_set_value(n_elts, 6, constant_val, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a tensor-valued quantity for a list of elements
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_tensor_by_val(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         dense_output,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *context,
                           cs_real_t                   *eval)
{
  CS_UNUSED(quant);
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  const cs_real_3_t  *constant_val = (const cs_real_3_t *)context;
  assert(eval != nullptr && constant_val != nullptr);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_tensor_on_subset(n_elts, elt_ids, constant_val, eval);
  else
    cs_array_real_set_tensor(n_elts, constant_val, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a scalar-valued quantity with only a time-dependent
 *        variation for a list of elements
 *        This function complies with the generic function type defined as
 *        cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_scalar_by_time_func(cs_lnum_t                   n_elts,
                                 const cs_lnum_t            *elt_ids,
                                 bool                        dense_output,
                                 const cs_mesh_t            *mesh,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant,
                                 cs_real_t                   time_eval,
                                 void                       *context,
                                 cs_real_t                  *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);

  cs_xdef_time_func_context_t  *tfc = (cs_xdef_time_func_context_t *)context;
  assert(tfc != nullptr);

  /* Evaluate the quantity only once */

  cs_real_t  _eval;
  tfc->func(time_eval, tfc->input, &_eval);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_scalar_on_subset(n_elts, elt_ids, _eval, eval);
  else
    cs_array_real_set_scalar(n_elts, _eval, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a vector-valued quantity with only a time-dependent
 *        variation for a list of elements.
 *        This function complies with the generic function type defined as
 *        cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_vector_by_time_func(cs_lnum_t                   n_elts,
                                 const cs_lnum_t            *elt_ids,
                                 bool                        dense_output,
                                 const cs_mesh_t            *mesh,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant,
                                 cs_real_t                   time_eval,
                                 void                       *context,
                                 cs_real_t                  *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);

  cs_xdef_time_func_context_t  *tfc = (cs_xdef_time_func_context_t *)context;
  assert(tfc != nullptr);

  /* Evaluate the quantity */

  cs_real_t  _eval[3];
  tfc->func(time_eval, tfc->input, _eval);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_vector_on_subset(n_elts, elt_ids, _eval, eval);
  else
    cs_array_real_set_vector(n_elts, _eval, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a tensor-valued quantity with a symmetric storage and with
 *        only a time-dependent variation for a list of elements
 *        This function complies with the generic function type defined as
 *        cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_symtens_by_time_func(cs_lnum_t                   n_elts,
                                  const cs_lnum_t            *elt_ids,
                                  bool                        dense_output,
                                  const cs_mesh_t            *mesh,
                                  const cs_cdo_connect_t     *connect,
                                  const cs_cdo_quantities_t  *quant,
                                  cs_real_t                   time_eval,
                                  void                       *context,
                                  cs_real_t                  *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);

  cs_xdef_time_func_context_t  *tfc = (cs_xdef_time_func_context_t *)context;
  assert(tfc != nullptr);

  /* Evaluate the quantity */

  cs_real_t  _eval[6];
  tfc->func(time_eval, tfc->input, _eval);

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_value_on_subset(n_elts, 6, elt_ids, _eval, eval);
  else
    cs_array_real_set_value(n_elts, 6, _eval, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a tensor-valued quantity with only a time-dependent
 *        variation for a list of elements
 *        This function complies with the generic function type defined as
 *        cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_tensor_by_time_func(cs_lnum_t                   n_elts,
                                 const cs_lnum_t            *elt_ids,
                                 bool                        dense_output,
                                 const cs_mesh_t            *mesh,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant,
                                 cs_real_t                   time_eval,
                                 void                       *context,
                                 cs_real_t                  *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);

  cs_xdef_time_func_context_t  *tfc = (cs_xdef_time_func_context_t *)context;
  assert(tfc != nullptr);

  /* Evaluate the quantity */

  cs_real_t  _eval[9];
  tfc->func(time_eval, tfc->input, _eval);
  cs_real_t   tens[3][3] = {{_eval[0], _eval[1], _eval[2]},
                            {_eval[3], _eval[4], _eval[5]},
                            {_eval[6], _eval[7], _eval[8]}};

  if (elt_ids != nullptr && !dense_output)
    cs_array_real_set_tensor_on_subset(n_elts, elt_ids, tens, eval);
  else
    cs_array_real_set_tensor(n_elts, tens, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at cells using an analytic function
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_cells_by_analytic(cs_lnum_t                    n_elts,
                                  const cs_lnum_t             *elt_ids,
                                  bool                         dense_output,
                                  const cs_mesh_t             *mesh,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  cs_real_t                    time_eval,
                                  void                        *context,
                                  cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  const cs_real_t *cell_centers
    = (quant != nullptr) ? quant->cell_centers : nullptr;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != nullptr);

  /* Evaluate the function for this time at cell centers */

  cx->func(time_eval, n_elts, elt_ids, cell_centers, dense_output, cx->input,
           eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at interior faces using an analytic
 *         function.
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_i_faces_by_analytic(cs_lnum_t                    n_elts,
                                    const cs_lnum_t             *elt_ids,
                                    bool                         dense_output,
                                    const cs_mesh_t             *mesh,
                                    const cs_cdo_connect_t      *connect,
                                    const cs_cdo_quantities_t   *quant,
                                    cs_real_t                    time_eval,
                                    void                        *context,
                                    cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  if (n_elts == 0)
    return;

  const cs_real_t *if_centers
    = (quant != nullptr) ? quant->i_face_center : nullptr;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != nullptr && eval != nullptr);

  /* Evaluate the function for this time at the center of interior faces */

  cx->func(time_eval, n_elts, elt_ids, if_centers, dense_output, cx->input,
           eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at border faces using an analytic
 *         function
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_b_faces_by_analytic(cs_lnum_t                    n_elts,
                                    const cs_lnum_t             *elt_ids,
                                    bool                         dense_output,
                                    const cs_mesh_t             *mesh,
                                    const cs_cdo_connect_t      *connect,
                                    const cs_cdo_quantities_t   *quant,
                                    cs_real_t                    time_eval,
                                    void                        *context,
                                    cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  if (n_elts == 0)
    return;

  const cs_real_t *bf_centers
    = (quant != nullptr) ? quant->b_face_center : nullptr;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != nullptr && eval != nullptr);

  /* Evaluate the function for this time at the center of border faces */

  cx->func(time_eval, n_elts, elt_ids, bf_centers, dense_output, cx->input,
           eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at vertices using an analytic function
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_vertices_by_analytic(cs_lnum_t                    n_elts,
                                     const cs_lnum_t             *elt_ids,
                                     bool                         dense_output,
                                     const cs_mesh_t             *mesh,
                                     const cs_cdo_connect_t      *connect,
                                     const cs_cdo_quantities_t   *quant,
                                     cs_real_t                    time_eval,
                                     void                        *context,
                                     cs_real_t                   *eval)
{
  CS_UNUSED(connect);

  if (n_elts == 0)
    return;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(eval != nullptr || cx != nullptr);

  const cs_real_t  *v_coords;
  if (quant != nullptr)
    v_coords = quant->vtx_coord;
  else if (mesh != nullptr)
    v_coords = mesh->vtx_coord;
  else {
    v_coords = nullptr; /* avoid a compilation warning */
    bft_error(__FILE__, __LINE__, 0, "%s: No vertex coordinates available.",
              __func__);
  }

  /* Evaluate the function for this time at vertices */

  cx->func(time_eval, n_elts, elt_ids, v_coords, dense_output, cx->input,
           eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at primal cells using a function
 *         associated to dof (dof = degrees of freedom).
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_cells_by_dof_func(cs_lnum_t                    n_elts,
                                  const cs_lnum_t             *elt_ids,
                                  bool                         dense_output,
                                  const cs_mesh_t             *mesh,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  cs_real_t                    time_eval,
                                  void                        *context,
                                  cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_eval);

  cs_xdef_dof_context_t  *cx = (cs_xdef_dof_context_t *)context;
  assert(cx != nullptr);

  /* Values of the function are defined at the cells */

  if (cs_flag_test(cx->dof_location, cs_flag_primal_cell))
    cx->func(n_elts, elt_ids, dense_output, cx->input, eval);
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid location.\n", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at vertices using a function
 *         associated to dof (dof = degrees of freedom).
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_vertices_by_dof_func(cs_lnum_t                    n_elts,
                                     const cs_lnum_t             *elt_ids,
                                     bool                         dense_output,
                                     const cs_mesh_t             *mesh,
                                     const cs_cdo_connect_t      *connect,
                                     const cs_cdo_quantities_t   *quant,
                                     cs_real_t                    time_eval,
                                     void                        *context,
                                     cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_eval);

  cs_xdef_dof_context_t  *cx = (cs_xdef_dof_context_t *)context;
  assert(cx != nullptr);

  /* Values of the function are defined at vertices */

  if (cs_flag_test(cx->dof_location, cs_flag_primal_vtx))
    cx->func(n_elts, elt_ids, dense_output, cx->input, eval);
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid location.\n", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at boundary faces using a function
 *         associated to dof (dof = degrees of freedom).
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_b_faces_by_dof_func(cs_lnum_t                    n_elts,
                                    const cs_lnum_t             *elt_ids,
                                    bool                         dense_output,
                                    const cs_mesh_t             *mesh,
                                    const cs_cdo_connect_t      *connect,
                                    const cs_cdo_quantities_t   *quant,
                                    cs_real_t                    time_eval,
                                    void                        *context,
                                    cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_eval);

  cs_xdef_dof_context_t  *cx = (cs_xdef_dof_context_t *)context;
  assert(cx != nullptr);

  /* Values of the function are defined at the boundary faces */

  if (cs_flag_test(cx->dof_location, cs_flag_boundary_face))
    cx->func(n_elts, elt_ids, dense_output, cx->input, eval);
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid location.\n", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a scalar-valued quantity at cells defined by an array.
 *         Array is assumed to be interlaced.
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_scalar_at_cells_by_array(cs_lnum_t                    n_elts,
                                      const cs_lnum_t             *elt_ids,
                                      bool                         dense_output,
                                      const cs_mesh_t             *mesh,
                                      const cs_cdo_connect_t      *connect,
                                      const cs_cdo_quantities_t   *quant,
                                      cs_real_t                    time_eval,
                                      void                        *context,
                                      cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;
  assert(eval != nullptr || cx != nullptr);
  assert(cx->stride == 1);

  if (cs_flag_test(cx->value_location, cs_flag_primal_cell)) {

    if (elt_ids != nullptr) {

      if (cx->full_length) {

        if (!dense_output)
          cs_array_real_copy_subset(
            n_elts, 1, elt_ids, CS_ARRAY_SUBSET_INOUT, cx->values, eval);
        else
          cs_array_real_copy_subset(
            n_elts, 1, elt_ids, CS_ARRAY_SUBSET_IN, cx->values, eval);
      }
      else { /* Not full length. We assume that n_elts and elt_ids are
                associated to the same zone as the current definition */

        if (!dense_output)
          cs_array_real_copy_subset(
            n_elts, 1, elt_ids, CS_ARRAY_SUBSET_OUT, cx->values, eval);
        else
          cs_array_real_copy(n_elts, cx->values, eval);
      }
    }
    else
      cs_array_real_copy(n_elts, cx->values, eval);
  }
  else if (cs_flag_test(cx->value_location, cs_flag_primal_vtx))
    cs_reco_scalar_v2c(
      n_elts, elt_ids, connect->c2v, quant, cx->values, dense_output, eval);

  else if (cs_flag_test(cx->value_location, cs_flag_dual_cell_byc))
    cs_reco_scalar_vbyc2c(
      n_elts, elt_ids, connect->c2v, quant, cx->values, dense_output, eval);

  else if (cs_flag_test(cx->value_location, cs_flag_primal_edge_byc))
    cs_reco_scalar_ebyc2c(
      n_elts, elt_ids, connect->c2e, quant, cx->values, dense_output, eval);

  else
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid support for the input array",
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a nd-valued quantity at cells defined by an array.
 *         Array is assumed to be interlaced.
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_nd_at_cells_by_array(cs_lnum_t                  n_elts,
                                  const cs_lnum_t           *elt_ids,
                                  bool                       dense_output,
                                  const cs_mesh_t           *mesh,
                                  const cs_cdo_connect_t    *connect,
                                  const cs_cdo_quantities_t *quant,
                                  cs_real_t                  time_eval,
                                  void                      *context,
                                  cs_real_t                 *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);
  CS_UNUSED(connect); /* Only in debug mode for a check */

  if (n_elts < 1)
    return;

  cs_xdef_array_context_t *cx = (cs_xdef_array_context_t *)context;
  assert(eval != nullptr || cx != nullptr);

  const int stride = cx->stride;

  if (cs_flag_test(cx->value_location, cs_flag_primal_cell)) {

    if (elt_ids != nullptr) {

      if (cx->full_length) {

        if (!dense_output)
          cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                    CS_ARRAY_SUBSET_INOUT,
                                    cx->values,
                                    eval);
        else
          cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                    CS_ARRAY_SUBSET_IN,
                                    cx->values,
                                    eval);

      }
      else { /* Not full length. We assume that n_elts and elt_ids are
                associated to the same zone as the current definition */

        if (!dense_output)
          cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                    CS_ARRAY_SUBSET_OUT,
                                    cx->values,
                                    eval);
        else
          cs_array_real_copy(stride*n_elts, cx->values, eval);
      }
    }
    else
      cs_array_real_copy(stride * n_elts, cx->values, eval);
  }
  else if (cs_flag_test(cx->value_location, cs_flag_dual_face_byc)) {

    assert(connect != nullptr && quant != nullptr);

    const cs_adjacency_t  *adj = cx->adjacency;
    assert(adj == connect->c2e);

    if (elt_ids != nullptr && !dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t c_id = elt_ids[i];
        cs_reco_dfbyc_at_cell_center(
          c_id, adj, quant, cx->values, eval + 3 * c_id);
      }
    }
    else if (elt_ids != nullptr && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_dfbyc_at_cell_center(
          elt_ids[i], adj, quant, cx->values, eval + 3 * i);
    }
    else {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_dfbyc_at_cell_center(i, adj, quant, cx->values, eval + 3*i);

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid case for the input array", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at vertices using an array
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_vertices_by_array(cs_lnum_t                    n_elts,
                                  const cs_lnum_t             *elt_ids,
                                  bool                         dense_output,
                                  const cs_mesh_t             *mesh,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  cs_real_t                    time_eval,
                                  void                        *context,
                                  cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;
  assert(eval != nullptr || cx != nullptr);

  const int  stride = cx->stride;

  if (cs_flag_test(cx->value_location, cs_flag_primal_vtx) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);

  if (cx->full_length) {

    if (elt_ids != nullptr && !dense_output)
      cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                CS_ARRAY_SUBSET_INOUT,
                                cx->values,
                                eval);

    else if (elt_ids != nullptr && dense_output)
      cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                CS_ARRAY_SUBSET_IN,
                                cx->values,
                                eval);

    else {

      assert(elt_ids == nullptr);
      cs_array_real_copy(n_elts*stride, cx->values, eval);

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: TODO. BC defined by an array on a subset.\n", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a quantity defined at boundary faces using an array
 *        This function complies with the generic function type defined as
 *        cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_b_faces_by_array(cs_lnum_t                    n_elts,
                                 const cs_lnum_t             *elt_ids,
                                 bool                         dense_output,
                                 const cs_mesh_t             *mesh,
                                 const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *quant,
                                 cs_real_t                    time_eval,
                                 void                        *context,
                                 cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;
  assert(eval != nullptr || cx != nullptr);

  const int  stride = cx->stride;

  if (cs_flag_test(cx->value_location, cs_flag_boundary_face) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);

  if (cx->full_length) {

    if (elt_ids != nullptr && !dense_output)
      cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                CS_ARRAY_SUBSET_INOUT,
                                cx->values,
                                eval);

    else if (elt_ids != nullptr && dense_output)
      cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                CS_ARRAY_SUBSET_IN,
                                cx->values,
                                eval);

    else {

      assert(elt_ids == nullptr);
      cs_array_real_copy(n_elts*stride, cx->values, eval);

    }

  }
  else {

    assert(elt_ids != nullptr);
    if (dense_output)
      cs_array_real_copy(n_elts*stride, cx->values, eval);
    else
      cs_array_real_copy_subset(n_elts, stride, elt_ids,
                                CS_ARRAY_SUBSET_OUT,
                                cx->values,
                                eval);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity inside a cell defined using a field
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cell_by_field(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         dense_output,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *context,
                           cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);

  if (n_elts < 1)
    return;

  cs_field_t  *field = (cs_field_t *)context;
  assert(eval != nullptr || field != nullptr);

  const cs_real_t  *values = field->val;

  switch (cs_mesh_location_get_type(field->location_id)) {

  case CS_MESH_LOCATION_CELLS:
    if (elt_ids != nullptr && !dense_output)
      cs_array_real_copy_subset(n_elts, field->dim, elt_ids,
                                CS_ARRAY_SUBSET_INOUT,
                                values,
                                eval);
    else
      cs_array_real_copy_subset(n_elts, field->dim, elt_ids,
                                CS_ARRAY_SUBSET_IN,
                                values,
                                eval);
    break;

  case CS_MESH_LOCATION_VERTICES: /* One operates a reconstruction at the cell
                                     centers */
    if (field->dim > 1)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Dimension %d not handled for field \"%s\".",
                __func__, field->dim, field->name);

    cs_reco_scalar_v2c(n_elts, elt_ids, connect->c2v, quant,
                       values,
                       dense_output,
                       eval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid case for the field \"%s\"",
              __func__, field->name);
    break;

  } /* Switch on mesh location */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at border faces using an analytic
 *         function
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       nullptr or pointer to a context structure
 * \param[in]      qtype         quadrature type
 * \param[in]      dim           dimension of the analytic function return
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_avg_at_b_faces_by_analytic(cs_lnum_t                    n_elts,
                                        const cs_lnum_t             *elt_ids,
                                        bool                    dense_output,
                                        const cs_mesh_t             *mesh,
                                        const cs_cdo_connect_t      *connect,
                                        const cs_cdo_quantities_t   *quant,
                                        cs_real_t                    time_eval,
                                        void                        *context,
                                        cs_quadrature_type_t         qtype,
                                        int                          dim,
                                        cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  assert(connect != nullptr && quant != nullptr);

  cs_xdef_analytic_context_t *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != nullptr);

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(dim, qtype);

  const cs_adjacency_t  *f2e = connect->f2e;
  const cs_adjacency_t  *e2v = connect->e2v;
  const cs_real_t  *xv = quant->vtx_coord;

  if (elt_ids == nullptr) { /* All boundary faces are selected */

#   pragma omp parallel for if (quant->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t bf_id = 0; bf_id < quant->n_b_faces; bf_id++) {

      const cs_lnum_t f_id = quant->n_i_faces + bf_id;
      const cs_quant_t pfq = cs_quant_set_face(f_id, quant);
      const cs_lnum_t  start = f2e->idx[f_id], end = f2e->idx[f_id+1];
      double *val_i = eval + dim*bf_id;

      /* Resetting */

      memset(val_i, 0, dim*sizeof(double));

      switch (end - start) {

      case CS_TRIANGLE_CASE:
        {
          cs_lnum_t v1, v2, v3;
          cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start,
                                         &v1, &v2, &v3);
          qfunc(time_eval, xv + 3*v1, xv + 3*v2, xv + 3*v3, pfq.meas,
                cx->func, cx->input, val_i);
        }
        break;

      default:
        for (cs_lnum_t j = start; j < end; j++) {

          const cs_lnum_t  _2e = 2*f2e->ids[j];
          const cs_lnum_t  v1 = e2v->ids[_2e];
          const cs_lnum_t  v2 = e2v->ids[_2e+1];

          qfunc(time_eval, xv + 3*v1, xv + 3*v2, pfq.center,
                cs_math_surftri(xv + 3*v1, xv + 3*v2, pfq.center),
                cx->func, cx->input, val_i);

        } /* Loop on edges */

      } /* Switch on the type of face. Special case for triangles */

      /* Compute the average */

      const double _os = 1./pfq.meas;
      for (int k = 0; k < dim; k++)
        val_i[k] *= _os;

    } /* Loop on faces */
  }
  else { /* There is an indirection list */

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) { /* Loop on selected faces */

      const cs_lnum_t  bf_id = elt_ids[i];
      const cs_lnum_t  f_id = quant->n_i_faces + bf_id;
      const cs_quant_t  pfq = cs_quant_set_face(f_id, quant);
      const cs_lnum_t  start = f2e->idx[f_id], end = f2e->idx[f_id+1];

      double  *val_i = dense_output ? eval + dim*i : eval + dim*bf_id;

      /* Resetting */

      memset(val_i, 0, dim*sizeof(double));

      switch (end - start) {

      case CS_TRIANGLE_CASE:
        {
          /* If triangle, one-shot computation */

          cs_lnum_t v1, v2, v3;
          cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start,
                                         &v1, &v2, &v3);
          qfunc(time_eval, xv + 3*v1, xv + 3*v2, xv + 3*v3,
                pfq.meas, cx->func, cx->input, val_i);
        }
        break;

      default:
        for (cs_lnum_t j = start; j < end; j++) {

          const cs_lnum_t  _2e = 2*f2e->ids[j];
          const cs_lnum_t  v1 = e2v->ids[_2e];
          const cs_lnum_t  v2 = e2v->ids[_2e+1];

          qfunc(time_eval, xv + 3*v1, xv + 3*v2, pfq.center,
                cs_math_surftri(xv + 3*v1, xv + 3*v2, pfq.center),
                cx->func, cx->input, val_i);

        } /* Loop on edges */

      } /* Switch on the type of face. Special case for triangles */

      /* Compute the average */

      const double _os = 1./pfq.meas;
      for (int k = 0; k < dim; k++)
        val_i[k] *= _os;

    } /* Loop on selected faces */

  }
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
