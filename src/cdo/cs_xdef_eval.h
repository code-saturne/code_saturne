#ifndef __CS_XDEF_EVAL_H__
#define __CS_XDEF_EVAL_H__

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
(cs_xdef_eval_cw_t) (const cs_cell_mesh_t    *cm,
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
(cs_xdef_eval_cw_xyz_t) (const cs_cell_mesh_t    *cm,
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
(cs_xdef_eval_cw_int_t) (const cs_cell_mesh_t    *cm,
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
(cs_xdef_eval_cw_face_t) (const cs_cell_mesh_t     *cm,
                          short int                 f,
                          cs_real_t                 time_eval,
                          void                     *input,
                          cs_quadrature_type_t      qtype,
                          cs_real_t                *eval);

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
 * \brief  Evaluate a scalar-valued quantity by a cellwise process
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_scalar_by_val(const cs_cell_mesh_t     *cm,
                              cs_real_t                 time_eval,
                              void                     *input,
                              cs_real_t                *eval);

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
 * \brief  Evaluate a vector-valued quantity by a cellwise process
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_vector_by_val(const cs_cell_mesh_t     *cm,
                              cs_real_t                 time_eval,
                              void                     *input,
                              cs_real_t                *eval);

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
 * \brief  Evaluate a tensor-valued quantity by a cellwise process
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_tensor_by_val(const cs_cell_mesh_t     *cm,
                              cs_real_t                 time_eval,
                              void                     *input,
                              cs_real_t                *eval);

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
                                        const short int              dim,
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
 * \brief  Evaluate a quantity defined using an analytic function by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_cell_by_analytic(const cs_cell_mesh_t       *cm,
                                 cs_real_t                   time_eval,
                                 void                       *input,
                                 cs_real_t                  *eval);

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
 * \brief  Evaluate a quantity at cells defined by an array.
 *         Array is assumed to be interlaced.
 *         Variation using a cs_cell_mesh_t structure
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_cell_by_array(const cs_cell_mesh_t      *cm,
                              cs_real_t                  time_eval,
                              void                      *input,
                              cs_real_t                 *eval);

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
/*!
 * \brief  Evaluate a quantity inside a cell defined using a field
 *         Variation using a cs_cell_mesh_t structure
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       value of the property at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_cell_by_field(const cs_cell_mesh_t        *cm,
                              cs_real_t                    time_eval,
                              void                        *input,
                              cs_real_t                   *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by analytic
 *         function at a precise location inside a cell
 *         Use of a cs_cell_mesh_t structure.
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  n_points   number of points where to compute the evaluation
 * \param[in]  xyz        where to compute the evaluation
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_at_xyz_by_analytic(const cs_cell_mesh_t       *cm,
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
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  n_points   number of points where to compute the evaluation
 * \param[in]  xyz        where to compute the evaluation
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_vector_at_xyz_by_val(const cs_cell_mesh_t       *cm,
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
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  n_points   number of points where to compute the evaluation
 * \param[in]  xyz        where to compute the evaluation
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_3_at_xyz_by_array(const cs_cell_mesh_t       *cm,
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
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  n_points   number of points where to compute the evaluation
 * \param[in]  xyz        where to compute the evaluation
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_3_at_xyz_by_field(const cs_cell_mesh_t    *cm,
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
cs_xdef_eval_cw_at_vtx_flux_by_val(const cs_cell_mesh_t     *cm,
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
cs_xdef_eval_cw_at_vtx_flux_by_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *input,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval);

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

void
cs_xdef_eval_cw_flux_by_val(const cs_cell_mesh_t     *cm,
                            short int                 f,
                            cs_real_t                 time_eval,
                            void                     *input,
                            cs_real_t                *eval);

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

void
cs_xdef_eval_cw_tensor_flux_by_val(const cs_cell_mesh_t     *cm,
                                   short int                 f,
                                   cs_real_t                 time_eval,
                                   void                     *input,
                                   cs_real_t                *eval);

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
cs_xdef_eval_cw_flux_by_analytic(const cs_cell_mesh_t      *cm,
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
cs_xdef_eval_cw_tensor_flux_by_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *input,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  f          local face id
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[in]  qtype      level of quadrature to use
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_face_avg_scalar_by_analytic(const cs_cell_mesh_t   *cm,
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
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  f          local face id
 * \param[in]  time_eval  physical time at which one evaluates the term
 * \param[in]  input      pointer to an input structure
 * \param[in]  qtype      level of quadrature to use
 * \param[out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_face_avg_scalar_by_analytic(const cs_cell_mesh_t   *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_face_avg_vector_by_analytic(const cs_cell_mesh_t    *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_face_avg_tensor_by_analytic(const cs_cell_mesh_t    *cm,
                                            short int                f,
                                            cs_real_t                t_eval,
                                            void                    *input,
                                            cs_quadrature_type_t     qtype,
                                            cs_real_t               *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (cs_xdef_t structure) by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  qtype    quadrature type
 * \param[in]  input    pointer to an input structure
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_avg_scalar_by_analytic(const cs_cell_mesh_t     *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  qtype    quadrature type
 * \param[in]  input    pointer to an input structure
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_avg_vector_by_analytic(const cs_cell_mesh_t     *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  qtype    quadrature type
 * \param[in]  input    pointer to an input structure
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_avg_tensor_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *input,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval);

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
 * \brief  Function pointer for evaluating a the reduction by averages of a
 *         analytic function by a cellwise process (usage of a
 *         cs_cell_mesh_t structure) which is hinged on integrals
 *         (faces first, then cell DoFs)
 *
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  qtype    quadrature type
 * \param[in]  input    pointer to an input structure
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_cw_vect_avg_reduction_by_analytic(const cs_cell_mesh_t     *cm,
                                               cs_real_t                 t_eval,
                                               void                     *input,
                                               cs_quadrature_type_t      qtype,
                                               cs_real_t                *eval);

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine to integrate an analytic function over a cell
 *
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   time at which the function is evaluated
 * \param[in]  ana      analytic function to integrate
 * \param[in]  input    pointer to an input structure
 * \param[in]  qfunc    quadrature function to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_int_on_cell(const cs_cell_mesh_t           *cm,
                         double                          t_eval,
                         cs_analytic_func_t             *ana,
                         void                           *input,
                         cs_quadrature_tetra_integral_t *qfunc,
                         cs_real_t                      *eval)
{
  switch (cm->type) {

  case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      qfunc(t_eval, cm->xv, cm->xv+3, cm->xv+6, cm->xv+9, cm->vol_c,
            ana, input, eval);
    }
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < cm->n_fc; ++f) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_onethird * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; // #vertices (=#edges)
      const short int *f2e_ids = cm->f2e_ids + start;

      assert(n_vf > 2);
      switch(n_vf){

      case CS_TRIANGLE_CASE: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          const double  *xv0 = cm->xv + 3*v0;
          const double  *xv1 = cm->xv + 3*v1;
          const double  *xv2 = cm->xv + 3*v2;

          qfunc(t_eval, xv0, xv1, xv2, cm->xc, hf_coef * pfq.meas,
                ana, input, eval);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            qfunc(t_eval, xv0, xv1, pfq.center, cm->xc, hf_coef * tef[e],
                  ana, input, eval);
          }
        }
        break;

      } /* End of switch */
    } /* End of loop on faces */

  }
  break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine to integrate an analytic function over a face
 *
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  t_eval   time at which the function is evaluated
 * \param[in]  f        local face id
 * \param[in]  ana      analytic function to integrate
 * \param[in]  input    pointer to an input structure
 * \param[in]  qfunc    quadrature function to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_int_on_face(const cs_cell_mesh_t          *cm,
                         double                         t_eval,
                         short int                      f,
                         cs_analytic_func_t            *ana,
                         void                          *input,
                         cs_quadrature_tria_integral_t *qfunc,
                         cs_real_t                     *eval)
{
  const cs_quant_t  pfq = cm->face[f];
  const int  start = cm->f2e_idx[f];
  const int  end = cm->f2e_idx[f+1];
  const short int n_vf = end - start; // #vertices (=#edges)
  const short int *f2e_ids = cm->f2e_ids + start;

  switch (n_vf) {
    case CS_TRIANGLE_CASE:
      {
        short int  v0, v1, v2;
        cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

        qfunc(t_eval, cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, pfq.meas,
              ana, input, eval);
      }
      break;
    default:
      {
        const double *tef = cm->tef + start;
        for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

          // Edge-related variables
          const short int e0  = f2e_ids[e];
          const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
          const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

          qfunc(t_eval, xv0, xv1, pfq.center, tef[e], ana, input, eval);
        }
      }

  } // Switch
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a scalar
 *         function defined through a descriptor (cs_xdef_t structure) by a
 *         cellwise process (usage of a cs_cell_mesh_t structure)
 *
 * \param[in]  cm      pointer to a cs_cell_mesh_t structure
 * \param[in]  f       local face id
 * \param[in]  t_eval  physical time at which one evaluates the term
 * \param[in]  input   pointer to an input structure
 * \param[in]  qtype   level of quadrature to use
 * \param[out] eval    result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_avg_scalar_by_value(const cs_cell_mesh_t     *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_avg_scalar_by_array(const cs_cell_mesh_t       *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_drhm_by_analytic(const cs_cell_mesh_t       *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_avg_vector_by_value(const cs_cell_mesh_t     *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_avg_vector_by_array(const cs_cell_mesh_t     *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_avg_tensor_by_value(const cs_cell_mesh_t     *cm,
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
 * \param[in]  cm       pointer to a cs_cell_mesh_t structure
 * \param[in]  f        local face id
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  input    pointer to an input structure
 * \param[in]  qtype    level of quadrature to use
 * \param[out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_xdef_eval_cw_face_avg_tensor_by_array(const cs_cell_mesh_t     *cm,
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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_EVAL_H__ */
