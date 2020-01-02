#ifndef __CS_XDEF_EVAL_H__
#define __CS_XDEF_EVAL_H__

/*============================================================================
 * Manage the (generic) evaluation of extended definitions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include "cs_cdo_quantities.h"
#include "cs_mesh.h"
#include "cs_quadrature.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Function pointer type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure)
 *
 * \param[in]      n_elts     number of elements to consider
 * \param[in]      elt_ids    list of element ids
 * \param[in]      compact    indirection for output (true or false)
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      input      pointer to an input structure cast on-the_fly
 * \param[in, out] eval       array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_xdef_eval_t) (cs_lnum_t                    n_elts,
                  const cs_lnum_t             *elt_ids,
                  bool                         compact,
                  const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  cs_real_t                    time_eval,
                  void                        *input,
                  cs_real_t                   *eval);

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
                           cs_real_t                   *eval);

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
                           cs_real_t                   *eval);

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
                           cs_real_t                   *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a scalar-valued quantity with only a time-dependent
 *        variation for a list of elements
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
cs_xdef_eval_scalar_at_cells_by_time_func(cs_lnum_t                   n_elts,
                                          const cs_lnum_t            *elt_ids,
                                          bool                        compact,
                                          const cs_mesh_t            *mesh,
                                          const cs_cdo_connect_t     *connect,
                                          const cs_cdo_quantities_t  *quant,
                                          cs_real_t                   time_eval,
                                          void                       *input,
                                          cs_real_t                  *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a vector-valued quantity with only a time-dependent
 *        variation for a list of elements
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
cs_xdef_eval_vector_at_cells_by_time_func(cs_lnum_t                   n_elts,
                                          const cs_lnum_t            *elt_ids,
                                          bool                        compact,
                                          const cs_mesh_t            *mesh,
                                          const cs_cdo_connect_t     *connect,
                                          const cs_cdo_quantities_t  *quant,
                                          cs_real_t                   time_eval,
                                          void                       *input,
                                          cs_real_t                  *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a tensor-valued quantity with only a time-dependent
 *        variation for a list of elements
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
cs_xdef_eval_tensor_at_cells_by_time_func(cs_lnum_t                   n_elts,
                                          const cs_lnum_t            *elt_ids,
                                          bool                        compact,
                                          const cs_mesh_t            *mesh,
                                          const cs_cdo_connect_t     *connect,
                                          const cs_cdo_quantities_t  *quant,
                                          cs_real_t                   time_eval,
                                          void                       *input,
                                          cs_real_t                  *eval);

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
                                  cs_real_t                   *eval);

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
                                    cs_real_t                   *eval);

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
                                        cs_real_t                   *eval);

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
                                     cs_real_t                   *eval);

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
                                      cs_real_t                   *eval);

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
                                  cs_real_t                   *eval);

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
                                  cs_real_t                   *eval);

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
                                        cs_real_t                  *eval);

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
                           cs_real_t                   *eval);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_EVAL_H__ */
