/*============================================================================
 * Manage the (generic) evaluation of extended definitions
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  const cs_real_t  *constant_val = (cs_real_t *)context;

  assert(eval != NULL && constant_val != NULL);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      eval[elt_ids[i]] = constant_val[0];

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      eval[i] = constant_val[0];

  }
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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  const cs_real_t  *constant_val = (cs_real_t *)context;

  assert(eval != NULL && constant_val != NULL);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *_res = eval + 3*elt_ids[i];

      _res[0] = constant_val[0];
      _res[1] = constant_val[1];
      _res[2] = constant_val[2];

    }

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *_res = eval + 3*i;

      _res[0] = constant_val[0];
      _res[1] = constant_val[1];
      _res[2] = constant_val[2];

    }

  }
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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  const cs_real_t  *constant_val = (const cs_real_t *)context;

  assert(eval != NULL && constant_val != NULL);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *shift_eval = eval + 6*elt_ids[i];
      for (int k = 0; k < 6; k++)
        shift_eval[k] = constant_val[k];

    }

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *shift_eval = eval + 6*i;
      for (int k = 0; k < 6; k++)
        shift_eval[k] = constant_val[k];

    }

  }
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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  const cs_real_3_t  *constant_val = (const cs_real_3_t *)context;

  assert(eval != NULL && constant_val != NULL);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *shift_eval = eval + 9*elt_ids[i];
      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          shift_eval[3*ki+kj] = constant_val[ki][kj];

    }

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *shift_eval = eval + 9*i;
      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          shift_eval[3*ki+kj] = constant_val[ki][kj];

    }

  }
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
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_scalar_at_cells_by_time_func(cs_lnum_t                   n_elts,
                                          const cs_lnum_t            *elt_ids,
                                          bool                   dense_output,
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
  assert(tfc != NULL);

  /* Evaluate the quantity only once */

  cs_real_t  _eval;
  tfc->func(time_eval, tfc->input, &_eval);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      eval[elt_ids[i]] = _eval;

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++)
      eval[i] = _eval;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a vector-valued quantity with only a time-dependent
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
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_vector_at_cells_by_time_func(cs_lnum_t                   n_elts,
                                          const cs_lnum_t            *elt_ids,
                                          bool                   dense_output,
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
  assert(tfc != NULL);

  /* Evaluate the quantity */

  cs_real_t  _eval[3];
  tfc->func(time_eval, tfc->input, _eval);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 3; k++)
        eval[3*elt_ids[i] + k] = _eval[k];

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 3; k++)
        eval[3*i+k] = _eval[k];

  }
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
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_symtens_at_cells_by_time_func(cs_lnum_t                  n_elts,
                                           const cs_lnum_t           *elt_ids,
                                           bool                   dense_output,
                                           const cs_mesh_t           *mesh,
                                           const cs_cdo_connect_t    *connect,
                                           const cs_cdo_quantities_t *quant,
                                           cs_real_t                  time_eval,
                                           void                      *context,
                                           cs_real_t                 *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);

  cs_xdef_time_func_context_t  *tfc = (cs_xdef_time_func_context_t *)context;
  assert(tfc != NULL);

  /* Evaluate the quantity */

  cs_real_t  _eval[6];
  tfc->func(time_eval, tfc->input, _eval);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 6; k++)
        eval[6*elt_ids[i] + k] = _eval[k];

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 6; k++)
        eval[6*i+k] = _eval[k];

  }
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
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_tensor_at_cells_by_time_func(cs_lnum_t                   n_elts,
                                          const cs_lnum_t            *elt_ids,
                                          bool                   dense_output,
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
  assert(tfc != NULL);

  /* Evaluate the quantity */

  cs_real_t  _eval[9];
  tfc->func(time_eval, tfc->input, _eval);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 9; k++)
        eval[9*elt_ids[i] + k] = _eval[k];

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 9; k++)
        eval[9*i+k] = _eval[k];

  }
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
 * \param[in]      context       NULL or pointer to a context structure
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

  const cs_real_t *cell_centers = (quant != NULL) ? quant->cell_centers : NULL;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != NULL);

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
 * \param[in]      context       NULL or pointer to a context structure
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

  const cs_real_t *if_centers = (quant != NULL) ? quant->i_face_center : NULL;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != NULL && eval != NULL);

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
 * \param[in]      context       NULL or pointer to a context structure
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

  const cs_real_t *bf_centers = (quant != NULL) ? quant->b_face_center : NULL;

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != NULL && eval != NULL);

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
 * \param[in]      context       NULL or pointer to a context structure
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
  assert(eval != NULL || cx != NULL);

  const cs_real_t  *v_coords;
  if (quant != NULL)
    v_coords = quant->vtx_coord;
  else if (mesh != NULL)
    v_coords = mesh->vtx_coord;
  else {
    v_coords = NULL;/* avoid a compilation warning */
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
 * \param[in]      context       NULL or pointer to a context structure
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
  assert(cx != NULL);

  /* Values of the function are defined at the cells */

  if (cs_flag_test(cx->loc, cs_flag_primal_cell))
    cx->func(n_elts, elt_ids, dense_output, cx->input,
             eval);
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
 * \param[in]      context       NULL or pointer to a context structure
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
  assert(cx != NULL);

  /* Values of the function are defined at vertices */

  if (cs_flag_test(cx->loc, cs_flag_primal_vtx))
    cx->func(n_elts, elt_ids, dense_output, cx->input,
             eval);
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
 * \param[in]      context       NULL or pointer to a context structure
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
  assert(cx != NULL);

  /* Values of the function are defined at the boundary faces */

  if (cs_flag_test(cx->loc, cs_flag_boundary_face))
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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;
  assert(eval != NULL || cx != NULL);
  assert(cx->stride == 1);

  if (cs_flag_test(cx->loc, cs_flag_primal_cell)) {

    if (elt_ids != NULL && !dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        eval[c_id] = cx->values[c_id];
      }

    }
    else if (elt_ids != NULL && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        eval[i] = cx->values[elt_ids[i]];

    }
    else {

      assert(elt_ids == NULL);
      memcpy(eval, cx->values, n_elts * sizeof(cs_real_t));

    }

  }
  else if (cs_flag_test(cx->loc, cs_flag_primal_vtx)) {

    assert(connect != NULL && quant != NULL);
    if (elt_ids != NULL && !dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_pv_at_cell_center(c_id, connect->c2v, quant, cx->values,
                                  eval + c_id);
      }

    }
    else if (elt_ids != NULL && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_pv_at_cell_center(elt_ids[i], connect->c2v, quant, cx->values,
                                  eval + i);

    }
    else {

      assert(elt_ids == NULL);
      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_pv_at_cell_center(i, connect->c2v, quant, cx->values,
                                  eval + i);

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
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
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_nd_at_cells_by_array(cs_lnum_t                    n_elts,
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

  if (n_elts == 0)
    return;

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;
  assert(eval != NULL || cx != NULL);

  const int  stride = cx->stride;

  if (cs_flag_test(cx->loc, cs_flag_primal_cell)) {

    assert(stride > 1);
    if (elt_ids != NULL && !dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < stride; k++)
          eval[stride*c_id + k] = cx->values[stride*c_id + k];
      }

    }
    else if (elt_ids != NULL && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < stride; k++)
          eval[stride*i + k] = cx->values[stride*c_id + k];
      }

    }
    else { /* All elements are considered */

      assert(elt_ids == NULL);
      memcpy(eval, cx->values, stride*n_elts * sizeof(cs_real_t));

    }

  }
  else if (cs_flag_test(cx->loc, cs_flag_dual_face_byc)) {

    assert(stride == 3);
    assert(connect!= NULL && quant != NULL);
    assert(cx->index == connect->c2e->idx);

    if (elt_ids != NULL && !dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_dfbyc_at_cell_center(c_id, connect->c2e, quant, cx->values,
                                     eval + c_id*stride);
      }

    }
    else if (elt_ids != NULL && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_dfbyc_at_cell_center(elt_ids[i], connect->c2e, quant,
                                     cx->values,
                                     eval + i*stride);

    }
    else {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_dfbyc_at_cell_center(i, connect->c2e, quant, cx->values,
                                     eval + i*stride);

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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;
  assert(eval != NULL || cx != NULL);

  const int  stride = cx->stride;

  if (cs_flag_test(cx->loc, cs_flag_primal_vtx)) {

    if (elt_ids != NULL && !dense_output) {

      switch (stride) {

      case 1: /* Scalar-valued */
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  v_id = elt_ids[i];
          eval[v_id] = cx->values[v_id];
        }
        break;

      default:
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  v_id = elt_ids[i];
          for (int j = 0; j < stride; j++)
            eval[stride*v_id + j] = cx->values[stride*v_id+j];
        }
        break;

      } /* End of switch */

    }
    else if (elt_ids != NULL && dense_output) {

      switch (stride) {

      case 1: /* Scalar-valued */
        for (cs_lnum_t i = 0; i < n_elts; i++)
          eval[i] = cx->values[elt_ids[i]];
        break;

      default:
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          for (int j = 0; j < stride; j++)
            eval[stride*i + j] = cx->values[stride*elt_ids[i] + j];
        }
        break;

      } /* End of switch */

    }
    else {

      assert(elt_ids == NULL);
      memcpy(eval, cx->values, n_elts*stride * sizeof(cs_real_t));

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
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
 * \param[in]      context       NULL or pointer to a context structure
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

  if (n_elts == 0)
    return;

  cs_field_t  *field = (cs_field_t *)context;
  assert(eval != NULL || field != NULL);

  cs_real_t  *values = field->val;

  const int  c_ml_id = cs_mesh_location_get_id_by_name(N_("cells"));
  const int  v_ml_id = cs_mesh_location_get_id_by_name(N_("vertices"));

  if (field->location_id == c_ml_id) {

    if (elt_ids != NULL && !dense_output) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < field->dim; k++)
          eval[field->dim*c_id + k] = values[field->dim*c_id + k];
      }
    }
    else if (elt_ids != NULL && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < field->dim; k++)
          eval[field->dim*i + k] = values[field->dim*c_id + k];
      }

    }
    else {

      assert(elt_ids == NULL);
      memcpy(eval, values, field->dim*n_elts * sizeof(cs_real_t));

    }

  }
  else if (field->location_id == v_ml_id) {

    assert(connect != NULL);
    if (field->dim > 1)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid field dimension.", __func__);

    if (elt_ids != NULL && !dense_output) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {

        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_pv_at_cell_center(c_id,
                                  connect->c2v,
                                  quant,
                                  values,
                                  eval + c_id);

      }
    }
    else if (elt_ids != NULL && dense_output) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {

        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_pv_at_cell_center(c_id,
                                  connect->c2v,
                                  quant,
                                  values,
                                  eval + i);

      }

    }
    else {

      assert(elt_ids == NULL);
      for (cs_lnum_t c_id = 0; c_id < n_elts; c_id++)
        cs_reco_pv_at_cell_center(c_id,
                                  connect->c2v,
                                  quant,
                                  values,
                                  eval + c_id);

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid case for the input field", __func__);
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
 * \param[in]      context       NULL or pointer to a context structure
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
  assert(connect != NULL && quant != NULL);

  cs_xdef_analytic_context_t *cx = (cs_xdef_analytic_context_t *)context;
  assert(cx != NULL);

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(dim, qtype);

  const cs_adjacency_t  *f2e = connect->f2e;
  const cs_adjacency_t  *e2v = connect->e2v;
  const cs_real_t  *xv = quant->vtx_coord;

  if (elt_ids == NULL) { /* All boundary faces are selected */

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
