#ifndef __CS_XDEF_CW_EVAL_H__
#define __CS_XDEF_CW_EVAL_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_mesh.h"
#include "cs_quadrature.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro definition (unset at the end of file)
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Function pointer type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure)
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_xdef_cw_eval_t) (const cs_cell_mesh_t    *cm,
                     cs_real_t                time_eval,
                     void                    *input,
                     cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity at several locations in a
 *         cell defined through a descriptor (cs_xdef_t structure).
 *         The algorithm may use a cs_cell_mesh_t structure.
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  n_points   number of points where to compute the evaluation
 * \param[in]  xyz        where to compute the evaluation
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_xdef_cw_eval_xyz_t) (const cs_cell_mesh_t    *cm,
                         cs_lnum_t                n_points,
                         const cs_real_t         *xyz,
                         cs_real_t                time_eval,
                         void                    *input,
                         cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  qtype      quadrature type
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_xdef_cw_eval_int_t) (const cs_cell_mesh_t    *cm,
                         cs_real_t                time_eval,
                         void                    *input,
                         cs_quadrature_type_t     qtype,
                         cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure)
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  f          local face id
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  qtype      quadrature type
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_xdef_cw_eval_face_t) (const cs_cell_mesh_t     *cm,
                          short int                 f,
                          cs_real_t                 time_eval,
                          void                     *input,
                          cs_quadrature_type_t      qtype,
                          cs_real_t                *eval);

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a scalar-valued quantity by a cellwise process
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_scalar_by_val(const cs_cell_mesh_t     *cm,
                              cs_real_t                 time_eval,
                              void                     *input,
                              cs_real_t                *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(time_eval);

  cs_real_t  *constant_val = (cs_real_t *)input;
  *eval = constant_val[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a vector-valued quantity by a cellwise process
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_vector_by_val(const cs_cell_mesh_t     *cm,
                              cs_real_t                 time_eval,
                              void                     *input,
                              cs_real_t                *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(time_eval);

  const cs_real_t  *constant_val = (cs_real_t *)input;

  eval[0] = constant_val[0];
  eval[1] = constant_val[1];
  eval[2] = constant_val[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a tensor-valued quantity by a cellwise process
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_tensor_by_val(const cs_cell_mesh_t     *cm,
                              cs_real_t                 time_eval,
                              void                     *input,
                              cs_real_t                *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(time_eval);

  const cs_real_3_t  *constant_val = (const cs_real_3_t *)input;
  for (int ki = 0; ki < 3; ki++)
    for (int kj = 0; kj < 3; kj++)
      eval[3*ki+kj] = constant_val[ki][kj];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by analytic
 *         function at a precise location inside a cell
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_vector_at_xyz_by_val(const cs_cell_mesh_t       *cm,
                                     cs_lnum_t                   n_points,
                                     const cs_real_t            *xyz,
                                     cs_real_t                   time_eval,
                                     void                       *input,
                                     cs_real_t                  *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(xyz);
  CS_UNUSED(time_eval);

  const cs_real_t  *constant_val = (cs_real_t *)input;

  for (int i = 0; i < n_points; i++) {
    eval[3*i    ] = constant_val[0];
    eval[3*i + 1] = constant_val[1];
    eval[2*i + 2] = constant_val[2];
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by values.
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation (set inside)
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_flux_by_val(const cs_cell_mesh_t     *cm,
                            short int                 f,
                            cs_real_t                 time_eval,
                            void                     *input,
                            cs_real_t                *eval)
{
  CS_UNUSED(time_eval);

  /* Sanity check */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ));

  const cs_real_t  *flux = (cs_real_t *)input;
  const cs_quant_t  fq = cm->face[f];

  eval[f] = fq.meas * _dp3(fq.unitv, flux);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by values.
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation (set inside)
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_tensor_flux_by_val(const cs_cell_mesh_t     *cm,
                                   short int                 f,
                                   cs_real_t                 time_eval,
                                   void                     *input,
                                   cs_real_t                *eval)
{
  CS_UNUSED(time_eval);

  cs_real_t  *flux = (cs_real_t *)input;
  const cs_quant_t  fq = cm->face[f];

  cs_math_33_3_product((const cs_real_t (*)[3])flux, fq.unitv, eval);
  for (int k = 0; k < 3; k++)
    eval[3*f+k] *= fq.meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      f       local face id
 * \param[in]      t_eval  physical time at which one evaluates the term
 * \param[in]      input   pointer to an input structure
 * \param[in]      qtype   level of quadrature to use
 * \param[in, out] eval    result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_scalar_face_avg_by_value(const cs_cell_mesh_t     *cm,
                                         short int                 f,
                                         cs_real_t                 t_eval,
                                         void                     *input,
                                         cs_quadrature_type_t      qtype,
                                         cs_real_t                *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(t_eval);
  CS_UNUSED(f);
  CS_UNUSED(qtype);

  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Array storing the evaluation should be allocated before"
              " the call to this function.", __func__);
  assert(input != NULL);

  eval[0] = ((const cs_real_t *)input)[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_scalar_face_avg_by_array(const cs_cell_mesh_t       *cm,
                                         short int                   f,
                                         cs_real_t                   t_eval,
                                         void                       *input,
                                         cs_quadrature_type_t        qtype,
                                         cs_real_t                  *eval)
{
  CS_UNUSED(t_eval);
  CS_UNUSED(qtype);

  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Array storing the evaluation should be allocated before"
              " the call to this function.", __func__);

  const cs_xdef_array_input_t *array_input =
    (const cs_xdef_array_input_t *)input;

  assert(input != NULL);
  assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

  eval[0] = array_input->values[cm->f_ids[f]];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating at the center of the face a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *         Since it's only an evaluation, the functions works for any dimension
 *         (supposed that the function is well defined)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_face_drhm_by_analytic(const cs_cell_mesh_t       *cm,
                                      short int                   f,
                                      cs_real_t                   t_eval,
                                      void                       *input,
                                      cs_quadrature_type_t        qtype,
                                      cs_real_t                  *eval)
{
  CS_UNUSED(qtype);
  cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)input;

  anai->func(t_eval, 1, NULL, cm->face[f].center, false, anai->input, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a vector
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_vector_face_avg_by_value(const cs_cell_mesh_t     *cm,
                                         short int                 f,
                                         cs_real_t                 t_eval,
                                         void                     *input,
                                         cs_quadrature_type_t      qtype,
                                         cs_real_t                *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(f);
  CS_UNUSED(t_eval);
  CS_UNUSED(qtype);

  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Array storing the evaluation should be allocated before"
              " the call to this function.", __func__);

  assert(input != NULL);

  memcpy(eval, (const cs_real_t *)input, 3*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a vector
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_vector_face_avg_by_array(const cs_cell_mesh_t     *cm,
                                         short int                 f,
                                         cs_real_t                 t_eval,
                                         void                     *input,
                                         cs_quadrature_type_t      qtype,
                                         cs_real_t                *eval)
{
  CS_UNUSED(t_eval);
  CS_UNUSED(qtype);

  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Array storing the evaluation should be allocated before"
              " the call to this function.", __func__);

  const cs_xdef_array_input_t *array_input =
    (const cs_xdef_array_input_t *)input;

  assert(input != NULL);
  assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

  memcpy(eval, array_input->values + 3*cm->f_ids[f], 3*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a tensor
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_tensor_face_avg_by_value(const cs_cell_mesh_t     *cm,
                                         short int                 f,
                                         cs_real_t                 t_eval,
                                         void                     *input,
                                         cs_quadrature_type_t      qtype,
                                         cs_real_t                *eval)
{
  CS_UNUSED(cm);
  CS_UNUSED(f);
  CS_UNUSED(t_eval);
  CS_UNUSED(qtype);

  assert(input != NULL);
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Array storing the evaluation should be allocated before"
              " the call to this function.", __func__);

  const cs_real_3_t  *constant_val = (const cs_real_3_t *)input;
  for (int ki = 0; ki < 3; ki++)
    for (int kj = 0; kj < 3; kj++)
      eval[3*ki+kj] = constant_val[ki][kj];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a tensor
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_cw_eval_tensor_face_avg_by_array(const cs_cell_mesh_t     *cm,
                                         short int                 f,
                                         cs_real_t                 t_eval,
                                         void                     *input,
                                         cs_quadrature_type_t      qtype,
                                         cs_real_t                *eval)
{
  CS_UNUSED(t_eval);
  CS_UNUSED(qtype);

  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Array storing the evaluation should be allocated before"
              " the call to this function.", __func__);

  const cs_xdef_array_input_t *array_input =
    (const cs_xdef_array_input_t *)input;

  assert(input != NULL);
  assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

  memcpy(eval, array_input->values + 9*cm->f_ids[f], 9*sizeof(cs_real_t));
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Integrate an analytic function over a face
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which the function is evaluated
 * \param[in]      ana      analytic function to integrate
 * \param[in]      input    pointer to an input structure
 * \param[in]      qfunc    quadrature function to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_face_int(const cs_cell_mesh_t            *cm,
                         double                           t_eval,
                         short int                        f,
                         cs_analytic_func_t              *ana,
                         void                            *input,
                         cs_quadrature_tria_integral_t   *qfunc,
                         cs_real_t                       *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_scalar_face_avg_by_analytic(const cs_cell_mesh_t   *cm,
                                            short int               f,
                                            cs_real_t               time_eval,
                                            void                   *input,
                                            cs_quadrature_type_t    qtype,
                                            cs_real_t              *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_face_avg_by_analytic(const cs_cell_mesh_t    *cm,
                                            short int                f,
                                            cs_real_t                t_eval,
                                            void                    *input,
                                            cs_quadrature_type_t     qtype,
                                            cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      input    pointer to an input structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_tensor_face_avg_by_analytic(const cs_cell_mesh_t    *cm,
                                            short int                f,
                                            cs_real_t                t_eval,
                                            void                    *input,
                                            cs_quadrature_type_t     qtype,
                                            cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Integrate an analytic function over a cell
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which the function is evaluated
 * \param[in]      ana      analytic function to integrate
 * \param[in]      input    pointer to an input structure
 * \param[in]      qfunc    quadrature function to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_cell_int_by_analytic(const cs_cell_mesh_t            *cm,
                                     double                           t_eval,
                                     cs_analytic_func_t              *ana,
                                     void                            *input,
                                     cs_quadrature_tetra_integral_t  *qfunc,
                                     cs_real_t                       *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      input    pointer to an input structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_scalar_avg_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *input,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      input    pointer to an input structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_avg_vector_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *input,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      input    pointer to an input structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_avg_tensor_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *input,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity defined using an analytic function by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_analytic(const cs_cell_mesh_t       *cm,
                            cs_real_t                   time_eval,
                            void                       *input,
                            cs_real_t                  *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity at cells defined by an array.
 *         Array is assumed to be interlaced.
 *         Variation using a cs_cell_mesh_t structure
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_array(const cs_cell_mesh_t      *cm,
                         cs_real_t                  time_eval,
                         void                      *input,
                         cs_real_t                 *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity inside a cell defined using a field
 *         Variation using a cs_cell_mesh_t structure
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       value of the property at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_field(const cs_cell_mesh_t        *cm,
                         cs_real_t                    time_eval,
                         void                        *input,
                         cs_real_t                   *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by analytic
 *         function at a precise location (x, y, z) inside a cell
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_at_xyz_by_analytic(const cs_cell_mesh_t       *cm,
                                   cs_lnum_t                   n_points,
                                   const cs_real_t            *xyz,
                                   cs_real_t                   time_eval,
                                   void                       *input,
                                   cs_real_t                  *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by analytic
 *         function at a precise location inside a cell
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_at_xyz_by_array(const cs_cell_mesh_t       *cm,
                                       cs_lnum_t                   n_points,
                                       const cs_real_t            *xyz,
                                       cs_real_t                   time_eval,
                                       void                       *input,
                                       cs_real_t                  *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by a field
 *         at a precise location inside a cell
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_at_xyz_by_field(const cs_cell_mesh_t    *cm,
                                       cs_lnum_t                n_points,
                                       const cs_real_t         *xyz,
                                       cs_real_t                time_eval,
                                       void                    *input,
                                       cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by values. The normal flux is then added to each portion of
 *         face related to a vertex.
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in, out] eval       result of the evaluation (updated inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_at_vtx_by_val(const cs_cell_mesh_t     *cm,
                                   short int                 f,
                                   cs_real_t                 time_eval,
                                   void                     *input,
                                   cs_real_t                *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by analytic function. The normal flux is then added to each
 *         portion of face related to a vertex.
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       result of the evaluation (updated inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_at_vtx_by_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *input,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by analytic function.
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       result of the evaluation (set inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_by_analytic(const cs_cell_mesh_t      *cm,
                                 short int                  f,
                                 cs_real_t                  time_eval,
                                 void                      *input,
                                 cs_quadrature_type_t       qtype,
                                 cs_real_t                 *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by analytic function.
 *         Case of vector-valued quantities.
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       result of the evaluation (set inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_tensor_flux_by_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *input,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine to integrate an analytic function over a cell and its faces
 *
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  ana      analytic function to integrate
 * \param[in]  input    pointer to an input structure
 * \param[in]  dim      dimension of the function
 * \param[in]  q_tet    quadrature function to use on tetrahedra
 * \param[in]  q_tri    quadrature function to use on triangles
 * \param[out] c_int    result of the evaluation on the cell
 * \param[out] f_int    result of the evaluation on the faces
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_int_on_cell_faces(const cs_cell_mesh_t            *cm,
                               cs_real_t                        t_eval,
                               cs_analytic_func_t              *ana,
                               void                            *input,
                               const short int                  dim,
                               cs_quadrature_tetra_integral_t  *q_tet,
                               cs_quadrature_tria_integral_t   *q_tri,
                               cs_real_t                       *c_int,
                               cs_real_t                       *f_int);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the reduction by averages of a
 *         analytic function by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *         (faces first, then cell).
 *         Vector-valued case
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      input    pointer to an input structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vect_avg_reduction_by_analytic(const cs_cell_mesh_t     *cm,
                                               cs_real_t                 t_eval,
                                               void                     *input,
                                               cs_quadrature_type_t      qtype,
                                               cs_real_t                *eval);

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS

#endif /* __CS_XDEF_CW_EVAL_H__ */
