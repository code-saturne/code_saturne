/*============================================================================
 * Manage the (generic) evaluation of extended definitions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

static const char _err_empty_array[] =
  " %s: Array storing the evaluation should be allocated before the call"
  " to this function.";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a scalar-valued quantity for a list of elements
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_scalar_by_val(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         compact,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *input,
                           cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);
  assert(eval != NULL);

  const cs_real_t  *constant_val = (cs_real_t *)input;

  if (elt_ids != NULL && !compact) {

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
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_vector_by_val(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         compact,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *input,
                           cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);
  assert(eval != NULL);

  const cs_real_t  *constant_val = (cs_real_t *)input;

  if (elt_ids != NULL && !compact) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  id = elt_ids[i];
      eval[3*id  ] = constant_val[0];
      eval[3*id+1] = constant_val[1];
      eval[3*id+2] = constant_val[2];
    }

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      eval[3*i  ] = constant_val[0];
      eval[3*i+1] = constant_val[1];
      eval[3*i+2] = constant_val[2];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a tensor-valued quantity for a list of elements
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_tensor_by_val(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         compact,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *input,
                           cs_real_t                   *eval)
{
  CS_UNUSED(quant);
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(time_eval);

  assert(eval != NULL);

  const cs_real_3_t  *constant_val = (const cs_real_3_t *)input;

  if (elt_ids != NULL && !compact) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i];
      cs_real_t  *shift_eval = eval + 9*id;
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
 * \brief  Evaluate a quantity defined at cells using an analytic function
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_cells_by_analytic(cs_lnum_t                    n_elts,
                                  const cs_lnum_t             *elt_ids,
                                  bool                         compact,
                                  const cs_mesh_t             *mesh,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  cs_real_t                    time_eval,
                                  void                        *input,
                                  cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)input;

  /* Evaluate the function for this time at the cell center */
  anai->func(time_eval,
             n_elts, elt_ids, quant->cell_centers,
             compact, /* Is output compacted ? */
             anai->input,
             eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at border faces using an analytic
 *         function
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_b_faces_by_analytic(cs_lnum_t                    n_elts,
                                    const cs_lnum_t             *elt_ids,
                                    bool                         compact,
                                    const cs_mesh_t             *mesh,
                                    const cs_cdo_connect_t      *connect,
                                    const cs_cdo_quantities_t   *quant,
                                    cs_real_t                    time_eval,
                                    void                        *input,
                                    cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)input;

  /* Evaluate the function for this time at the border face center */
  anai->func(time_eval,
             n_elts, elt_ids, quant->b_face_center,
             compact,  // compacted output ?
             anai->input,
             eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined at border faces using an analytic
 *         function
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[in]  qtype      quadrature type
 * \param[in]  dim        dimension of the analytic function return
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_avg_at_b_faces_by_analytic(cs_lnum_t                    n_elts,
                                        const cs_lnum_t             *elt_ids,
                                        bool                         compact,
                                        const cs_mesh_t             *mesh,
                                        const cs_cdo_connect_t      *connect,
                                        const cs_cdo_quantities_t   *quant,
                                        cs_real_t                    time_eval,
                                        void                        *input,
                                        cs_quadrature_type_t         qtype,
                                        const int                    dim,
                                        cs_real_t                   *eval)
{
  CS_UNUSED(mesh);

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(dim, qtype);
  cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)input;

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

      switch (end - start) {

      case CS_TRIANGLE_CASE:
        {
          cs_lnum_t v1, v2, v3;
          cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start,
                                         &v1, &v2, &v3);
          qfunc(time_eval, xv + 3*v1, xv + 3*v2, xv + 3*v3, pfq.meas,
                anai->func, anai->input, val_i);
        }
        break;

      default:
        for (cs_lnum_t j = start; j < end; j++) {

          const cs_lnum_t  _2e = 2*f2e->ids[j];
          const cs_lnum_t  v1 = e2v->ids[_2e];
          const cs_lnum_t  v2 = e2v->ids[_2e+1];

          qfunc(time_eval, xv + 3*v1, xv + 3*v2, pfq.center,
                cs_math_surftri(xv + 3*v1, xv + 3*v2, pfq.center),
                anai->func, anai->input, val_i);

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

      double  *val_i = compact ? eval + dim*i : eval + dim*bf_id;

      switch (end - start) {

      case CS_TRIANGLE_CASE:
        {
          /* If triangle, one-shot computation */
          cs_lnum_t v1, v2, v3;
          cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start,
                                         &v1, &v2, &v3);
          qfunc(time_eval, xv + 3*v1, xv + 3*v2, xv + 3*v3,
                pfq.meas, anai->func, anai->input, val_i);
        }
        break;

      default:
        for (cs_lnum_t j = start; j < end; j++) {

          const cs_lnum_t  _2e = 2*f2e->ids[j];
          const cs_lnum_t  v1 = e2v->ids[_2e];
          const cs_lnum_t  v2 = e2v->ids[_2e+1];

          qfunc(time_eval, xv + 3*v1, xv + 3*v2, pfq.center,
                cs_math_surftri(xv + 3*v1, xv + 3*v2, pfq.center),
                anai->func, anai->input, val_i);

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
/*!
 * \brief  Evaluate a quantity defined at vertices using an analytic function
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_vertices_by_analytic(cs_lnum_t                    n_elts,
                                     const cs_lnum_t             *elt_ids,
                                     bool                         compact,
                                     const cs_mesh_t             *mesh,
                                     const cs_cdo_connect_t      *connect,
                                     const cs_cdo_quantities_t   *quant,
                                     cs_real_t                    time_eval,
                                     void                        *input,
                                     cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)input;

  /* Evaluate the function for this time at the cell center */
  anai->func(time_eval,
             n_elts, elt_ids, quant->vtx_coord,
             compact,  // compacted output ?
             anai->input,
             eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a scalar-valued quantity at cells defined by an array.
 *         Array is assumed to be interlaced.
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_scalar_at_cells_by_array(cs_lnum_t                    n_elts,
                                      const cs_lnum_t             *elt_ids,
                                      bool                         compact,
                                      const cs_mesh_t             *mesh,
                                      const cs_cdo_connect_t      *connect,
                                      const cs_cdo_quantities_t   *quant,
                                      cs_real_t                    time_eval,
                                      void                        *input,
                                      cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);

  cs_xdef_array_input_t  *array_input = (cs_xdef_array_input_t *)input;

  assert(array_input->stride == 1);

  if ((array_input->loc & cs_flag_primal_cell) == cs_flag_primal_cell) {

    if (elt_ids != NULL && !compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        eval[c_id] = array_input->values[c_id];
      }

    }
    else if (elt_ids != NULL && compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        eval[i] = array_input->values[elt_ids[i]];

    }
    else {

      assert(elt_ids == NULL);
      memcpy(eval, array_input->values, n_elts * sizeof(cs_real_t));

    }

  }
  else if ((array_input->loc & cs_flag_primal_vtx) == cs_flag_primal_vtx) {

    if (elt_ids != NULL && !compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_pv_at_cell_center(c_id,
                                  connect->c2v,
                                  quant,
                                  array_input->values,
                                  eval + c_id);
      }

    }
    else if (elt_ids != NULL && compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_pv_at_cell_center(elt_ids[i],
                                  connect->c2v,
                                  quant,
                                  array_input->values,
                                  eval + i);

    }
    else {

      assert(elt_ids == NULL);
      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_pv_at_cell_center(i,
                                  connect->c2v,
                                  quant,
                                  array_input->values,
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
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_nd_at_cells_by_array(cs_lnum_t                    n_elts,
                                  const cs_lnum_t             *elt_ids,
                                  bool                         compact,
                                  const cs_mesh_t             *mesh,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  cs_real_t                    time_eval,
                                  void                        *input,
                                  cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);

  cs_xdef_array_input_t  *array_input = (cs_xdef_array_input_t *)input;

  const int  stride = array_input->stride;

  if (cs_flag_test(array_input->loc, cs_flag_primal_cell)) {

    assert(stride > 1);
    if (elt_ids != NULL && !compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < stride; k++)
          eval[stride*c_id + k] = array_input->values[stride*c_id + k];
      }

    }
    else if (elt_ids != NULL && compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < stride; k++)
          eval[stride*i + k] = array_input->values[stride*c_id + k];
      }

    }
    else {

      assert(elt_ids == NULL);
      memcpy(eval, array_input->values, stride*n_elts * sizeof(cs_real_t));

    }

  }
  else if (cs_flag_test(array_input->loc, cs_flag_dual_face_byc)) {

    assert(stride == 3);
    assert(array_input->index == connect->c2e->idx);

    if (elt_ids != NULL && !compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_dfbyc_at_cell_center(c_id,
                                     connect->c2e,
                                     quant,
                                     array_input->values,
                                     eval + c_id*stride);
      }

    }
    else if (elt_ids != NULL && compact) {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_dfbyc_at_cell_center(elt_ids[i],
                                     connect->c2e,
                                     quant,
                                     array_input->values,
                                     eval + i*stride);

    }
    else {

      for (cs_lnum_t i = 0; i < n_elts; i++)
        cs_reco_dfbyc_at_cell_center(i,
                                     connect->c2e,
                                     quant,
                                     array_input->values,
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
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_vertices_by_array(cs_lnum_t                    n_elts,
                                  const cs_lnum_t             *elt_ids,
                                  bool                         compact,
                                  const cs_mesh_t             *mesh,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  cs_real_t                    time_eval,
                                  void                        *input,
                                  cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_eval);

  cs_xdef_array_input_t  *array_input = (cs_xdef_array_input_t *)input;

  const int  stride = array_input->stride;

  if (cs_flag_test(array_input->loc, cs_flag_primal_vtx)) {

    if (elt_ids != NULL && !compact) {

      switch (stride) {

      case 1: /* Scalar-valued */
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  v_id = elt_ids[i];
          eval[v_id] = array_input->values[v_id];
        }
        break;

      default:
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  v_id = elt_ids[i];
          for (int j = 0; j < stride; j++)
            eval[stride*v_id + j] = array_input->values[stride*v_id+j];
        }
        break;

      } /* End of switch */

    }
    else if (elt_ids != NULL && compact) {

      switch (stride) {

      case 1: /* Scalar-valued */
        for (cs_lnum_t i = 0; i < n_elts; i++)
          eval[i] = array_input->values[elt_ids[i]];
        break;

      default:
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          for (int j = 0; j < stride; j++)
            eval[stride*i + j] = array_input->values[stride*elt_ids[i] + j];
        }
        break;

      } /* End of switch */

    }
    else {

      assert(elt_ids == NULL);
      memcpy(eval, array_input->values, n_elts*stride * sizeof(cs_real_t));

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a vector-valued quantity at all vertices defined by an
 *         array.
 *         Array is assumed to be interlaced.
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_3_at_all_vertices_by_array(cs_lnum_t                   n_elts,
                                        const cs_lnum_t            *elt_ids,
                                        bool                        compact,
                                        const cs_mesh_t            *mesh,
                                        const cs_cdo_connect_t     *connect,
                                        const cs_cdo_quantities_t  *quant,
                                        cs_real_t                   time_eval,
                                        void                       *input,
                                        cs_real_t                  *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);
  CS_UNUSED(compact);

  cs_xdef_array_input_t  *array_input = (cs_xdef_array_input_t *)input;

  const int  stride = array_input->stride;

  if (elt_ids != NULL || n_elts < quant->n_vertices)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case\n", __func__);

  double  *dc_vol = NULL;
  BFT_MALLOC(dc_vol, quant->n_vertices, double);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
    dc_vol[i] = 0;

  if (cs_flag_test(array_input->loc, cs_flag_primal_cell)) {

    assert(stride == 3);
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Retrieve the cell vector */
      cs_real_3_t  cell_vector;
      for (int k = 0; k < stride; k++)
        cell_vector[k] = array_input->values[stride*c_id + k];

      /* Interpolate with a weighting related to the vertex volume in each
         cell */
      const cs_lnum_t  *c2v_idx = connect->c2v->idx + c_id;
      const cs_lnum_t  *c2v_ids = connect->c2v->ids + c2v_idx[0];
      const double  *vol_vc = quant->dcell_vol + c2v_idx[0];

      for (short int v = 0; v < c2v_idx[1]-c2v_idx[0]; v++) {

        const cs_lnum_t  v_id = c2v_ids[v];

        dc_vol[v_id] += vol_vc[v];
        cs_real_t  *v_val = eval + 3*v_id;
        for (int k = 0; k < 3; k++) v_val[k] += vol_vc[v] * cell_vector[k];

      } // Loop on cell vertices

    } // Loop on cells

#   pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

      const double  inv_dcvol = 1/dc_vol[v_id];
      cs_real_t *v_val = eval + 3*v_id;
      for (int k = 0; k < 3; k++) v_val[k] *= inv_dcvol;

    } // Loop on vertices

  }
  else if (cs_flag_test(array_input->loc, cs_flag_dual_face_byc)) {

    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Compute a estimated cell vector */
      cs_real_3_t  cell_vector;
      cs_reco_dfbyc_at_cell_center(c_id,
                                   connect->c2e,
                                   quant,
                                   array_input->values,
                                   cell_vector);

      /* Interpolate with a weighting related to the vertex volume in each
         cell */
      const cs_lnum_t  *c2v_idx = connect->c2v->idx + c_id;
      const cs_lnum_t  *c2v_ids = connect->c2v->ids + c2v_idx[0];
      const double  *vol_vc = quant->dcell_vol + c2v_idx[0];

      for (short int v = 0; v < c2v_idx[1]-c2v_idx[0]; v++) {

        const cs_lnum_t  v_id = c2v_ids[v];

        dc_vol[v_id] += vol_vc[v];
        cs_real_t  *v_val = eval + 3*v_id;
        for (int k = 0; k < 3; k++) v_val[k] += vol_vc[v] * cell_vector[k];

      } // Loop on cell vertices

    } // Loop on cells

#   pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

      const double  inv_dcvol = 1/dc_vol[v_id];
      cs_real_t *v_val = eval + 3*v_id;
      for (int k = 0; k < 3; k++) v_val[k] *= inv_dcvol;

    } // Loop on vertices

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid case for the input array", __func__);

  /* Free temporary buffer */
  BFT_FREE(dc_vol);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity inside a cell defined using a field
 *
 * \param[in]  n_elts     number of elements to consider
 * \param[in]  elt_ids    list of element ids
 * \param[in]  compact    true:no indirection, false:indirection for output
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cell_by_field(cs_lnum_t                    n_elts,
                           const cs_lnum_t             *elt_ids,
                           bool                         compact,
                           const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           cs_real_t                    time_eval,
                           void                        *input,
                           cs_real_t                   *eval)
{
  CS_UNUSED(mesh);
  CS_UNUSED(time_eval);

  cs_field_t  *field = (cs_field_t *)input;
  assert(field != NULL);
  cs_real_t  *values = field->val;
  const int  c_ml_id = cs_mesh_location_get_id_by_name(N_("cells"));
  const int  v_ml_id = cs_mesh_location_get_id_by_name(N_("vertices"));

  if (field->location_id == c_ml_id) {

    if (elt_ids != NULL && !compact) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        const cs_lnum_t  c_id = elt_ids[i];
        for (int k = 0; k < field->dim; k++)
          eval[field->dim*c_id + k] = values[field->dim*c_id + k];
      }
    }
    else if (elt_ids != NULL && compact) {

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

    if (field->dim > 1)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid field dimension.", __func__);

    if (elt_ids != NULL && !compact) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {

        const cs_lnum_t  c_id = elt_ids[i];
        cs_reco_pv_at_cell_center(c_id,
                                  connect->c2v,
                                  quant,
                                  values,
                                  eval + c_id);

      }
    }
    else if (elt_ids != NULL && compact) {

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

#undef _dp3

END_C_DECLS
