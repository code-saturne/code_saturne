#ifndef __CS_EVALUATE_H__
#define __CS_EVALUATE_H__

/*============================================================================
 * Functions and structures to deal with evaluation of quantities
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
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
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute reduced quantities for an array of size equal to dim * n_x
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in]      n_x     number of elements
 * \param[in]      array   array to analyze
 * \param[in]      w_x     weight to apply (may be set to  NULL)
 * \param[in, out] min     resulting min array (size: dim, or 4 if dim = 3)
 * \param[in, out] max     resulting max array (size: dim, or 4 if dim = 3)
 * \param[in, out] wsum    (weighted) sum array (size: dim, or 4 if dim = 3)
 * \param[in, out] asum    (weighted) sum of absolute values (same size as wsum)
 * \param[in, out] ssum    (weighted) sum of squared values (same size as wsum)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_array_reduction(int                     dim,
                            cs_lnum_t               n_x,
                            const cs_real_t        *array,
                            const cs_real_t        *w_x,
                            cs_real_t              *min,
                            cs_real_t              *max,
                            cs_real_t              *wsum,
                            cs_real_t              *asum,
                            cs_real_t              *ssum);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute reduced quantities for an array attached to either vertex,
 *         face or edge DoFs
 *         The weight to apply to each entity x is scanned using the adjacency
 *         structure. array size is equal to dim * n_x
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in]      n_x     number of elements
 * \param[in]      array   array to analyze
 * \param[in]      w_x     weight to apply (may be set to  NULL)
 * \param[in, out] min     resulting min array (size: dim, or 4 if dim = 3)
 * \param[in, out] max     resulting max array (size: dim, or 4 if dim = 3)
 * \param[in, out] vsum    (weighted) sum array (size: dim, or 4 if dim = 3)
 * \param[in, out] asum    (weighted) sum of absolute values (same size as vsum)
 * \param[in, out] ssum    (weighted) sum of squared values (same size as vsum)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_scatter_array_reduction(int                     dim,
                                    cs_lnum_t               n_x,
                                    const cs_real_t        *array,
                                    const cs_adjacency_t   *c2x,
                                    const cs_real_t        *w_x,
                                    cs_real_t              *min,
                                    cs_real_t              *max,
                                    cs_real_t              *wsum,
                                    cs_real_t              *asum,
                                    cs_real_t              *ssum);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted L2-norm of an array. The weight is scanned
 *         by a \ref cs_adjacency_t structure
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 * \param[in]  c2x     ajacency structure from cell to x entities (mandatory)
 * \param[in]  w_c2x   weight to apply (mandatory), scanned by c2x
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_square_wc2x_norm(const cs_real_t        *array,
                             const cs_adjacency_t   *c2x,
                             const cs_real_t        *w_c2x);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted L2-norm of an array. The weight is scanned
 *         by a \ref cs_adjacency_t structure.
 *         Case of a vector-valued array.
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 * \param[in]  c2x     ajacency structure from cell to x entities (mandatory)
 * \param[in]  w_c2x   weight to apply (mandatory), scanned by c2x
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_3_square_wc2x_norm(const cs_real_t        *array,
                               const cs_adjacency_t   *c2x,
                               const cs_real_t        *w_c2x);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the relative norm of the difference of two arrays scanned
 *         by the same \ref cs_adjacency_t structure. Normalization is done
 *         with the reference array.
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 * \param[in]  ref     array used for normalization and difference
 * \param[in]  c2x     ajacency structure from cell to x entities (mandatory)
 * \param[in]  w_c2x   weight to apply (mandatory), scanned by c2x
 *
 * \return the computed square weighted and normalized L2-norm of the
 *          difference between array and reference
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_delta_square_wc2x_norm(const cs_real_t        *array,
                                   const cs_real_t        *ref,
                                   const cs_adjacency_t   *c2x,
                                   const cs_real_t        *w_c2x);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the relative norm of the difference of two arrays scanned
 *         by the same \ref cs_adjacency_t structure. Normalization is done
 *         with the reference array.
 *         The quantities computed are synchronized in parallel.
 *         Case of vector-valued arrays.
 *
 * \param[in]  array   array to analyze
 * \param[in]  ref     array used for normalization and difference
 * \param[in]  c2x     ajacency structure from cell to x entities (mandatory)
 * \param[in]  w_c2x   weight to apply (mandatory), scanned by c2x
 *
 * \return the computed square weighted and normalized L2-norm of the
 *          difference between array and reference
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_delta_3_square_wc2x_norm(const cs_real_t        *array,
                                     const cs_real_t        *ref,
                                     const cs_adjacency_t   *c2x,
                                     const cs_real_t        *w_c2x);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         The value defined by the analytic function is by unity of volume
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      def         pointer to a cs_xdef_t structure
 * \param[in]      time_eval   physical time at which one evaluates the term
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_by_analytic(cs_flag_t           dof_flag,
                                const cs_xdef_t    *def,
                                cs_real_t           time_eval,
                                cs_real_t           retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity defined by a value in the case of a density
 *         field for all the degrees of freedom
 *         Accessor to the value is by unit of volume
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      def       pointer to a cs_xdef_t structure
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_by_value(cs_flag_t          dof_flag,
                             const cs_xdef_t   *def,
                             cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at vertices
 *         when the definition relies on an analytic expression
 *
 * \param[in]      def           pointer to a cs_xdef_t pointer
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      n_v_selected  number of selected vertices
 * \param[in]      selected_lst  list of selected vertices
 * \param[in, out] retval        pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_vertices_by_analytic(const cs_xdef_t   *def,
                                              const cs_real_t    time_eval,
                                              const cs_lnum_t    n_v_selected,
                                              const cs_lnum_t   *selected_lst,
                                              cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at face centers
 *         when the definition relies on an analytic expression
 *
 * \param[in]      def           pointer to a cs_xdef_t pointer
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      n_f_selected  number of selected faces
 * \param[in]      selected_lst  list of selected faces
 * \param[in, out] retval        pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_faces_by_analytic(const cs_xdef_t   *def,
                                           const cs_real_t    time_eval,
                                           const cs_lnum_t    n_f_selected,
                                           const cs_lnum_t   *selected_lst,
                                           cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at cell centers
 *         when the definition relies on an analytic expression
 *
 * \param[in]      def         pointer to a cs_xdef_t pointer
 * \param[in]      time_eval   physical time at which one evaluates the term
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_cells_by_analytic(const cs_xdef_t    *def,
                                           const cs_real_t     time_eval,
                                           cs_real_t           retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF in the case of a potential field in order
 *         to put a given quantity inside the volume associated to the zone
 *         related to the given definition
 *         wvals may be NULL.
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] vvals     pointer to the first array of computed values
 * \param[in, out] wvals     pointer to the second array of computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_by_qov(cs_flag_t          dof_flag,
                             const cs_xdef_t   *def,
                             cs_real_t          vvals[],
                             cs_real_t          wvals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field at vertices from a definition by a
 *         constant value
 *
 * \param[in]      def             pointer to a cs_xdef_t pointer
 * \param[in]      n_v_selected    number of selected vertices
 * \param[in]      selected_lst    list of selected vertices
 * \param[in, out] retval          pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_vertices_by_value(const cs_xdef_t   *def,
                                           const cs_lnum_t    n_v_selected,
                                           const cs_lnum_t   *selected_lst,
                                           cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field atface centers from a definition by a
 *         constant value
 *
 * \param[in]      def             pointer to a cs_xdef_t pointer
 * \param[in]      n_f_selected    number of selected faces
 * \param[in]      selected_lst    list of selected faces
 * \param[in, out] retval          pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_faces_by_value(const cs_xdef_t   *def,
                                        const cs_lnum_t    n_f_selected,
                                        const cs_lnum_t   *selected_lst,
                                        cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field at cell centers from a definition by
 *         value
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_cells_by_value(const cs_xdef_t   *def,
                                        cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the circulation along a selection of (primal) edges.
 *         Circulation is defined thanks to a constant vector field (by value)
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      n_e_selected   number of selected edges
 * \param[in]      selected_lst   list of selected edges
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_circulation_along_edges_by_value(const cs_xdef_t   *def,
                                             const cs_lnum_t    n_e_selected,
                                             const cs_lnum_t   *selected_lst,
                                             cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the circulation along a selection of (primal) edges.
 *         Circulation is defined thanks to an array
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      n_e_selected   number of selected edges
 * \param[in]      selected_lst   list of selected edges
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_circulation_along_edges_by_array(const cs_xdef_t   *def,
                                             const cs_lnum_t    n_e_selected,
                                             const cs_lnum_t   *selected_lst,
                                             cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the circulation along a selection of (primal) edges.
 *         Circulation is defined by an analytical function.
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      n_e_selected   number of selected edges
 * \param[in]      selected_lst   list of selected edges
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_circulation_along_edges_by_analytic(const cs_xdef_t   *def,
                                                const cs_real_t    time_eval,
                                                const cs_lnum_t    n_e_selected,
                                                const cs_lnum_t   *selected_lst,
                                                cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the faces
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      n_f_selected   number of selected faces
 * \param[in]      selected_lst   list of selected faces
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_faces_by_value(const cs_xdef_t   *def,
                                      const cs_lnum_t    n_f_selected,
                                      const cs_lnum_t   *selected_lst,
                                      cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the faces.
 *         Warning: retval has to be initialize before calling this function.
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      n_f_selected   number of selected faces
 * \param[in]      selected_lst   list of selected faces
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_faces_by_analytic(const cs_xdef_t   *def,
                                         const cs_real_t    time_eval,
                                         const cs_lnum_t    n_f_selected,
                                         const cs_lnum_t   *selected_lst,
                                         cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells_by_value(const cs_xdef_t   *def,
                                      cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells_by_array(const cs_xdef_t   *def,
                                      cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells
 *         Warning: retval has to be initialize before calling this function.
 *
 * \param[in]      def        pointer to a cs_xdef_t pointer
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] retval     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells_by_analytic(const cs_xdef_t   *def,
                                         cs_real_t          time_eval,
                                         cs_real_t          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the integral over the full computational domain of a
 *         quantity defined by an array
 *
 * \param[in]      array_loc  flag indicating where are located values
 * \param[in]      array_val  array of values
 *
 * \return the value of the integration
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_scal_domain_integral_by_array(cs_flag_t         array_loc,
                                          const cs_real_t  *array_val);

/*============================================================================
 * Inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the faces
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      n_f_selected   number of selected faces
 * \param[in]      selected_lst   list of selected faces
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_evaluate_average_on_faces(const cs_xdef_t   *def,
                             cs_real_t          time_eval,
                             const cs_lnum_t    n_f_selected,
                             const cs_lnum_t   *selected_lst,
                             cs_real_t          retval[])
{
  /* Sanity checks */
  assert(def != NULL);

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    cs_evaluate_average_on_faces_by_value(def,
                                          n_f_selected, selected_lst,
                                          retval);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_average_on_faces_by_analytic(def,
                                             time_eval,
                                             n_f_selected, selected_lst,
                                             retval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Case not handled yet.", __func__);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells
 *
 * \param[in]      def        pointer to a cs_xdef_t pointer
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] retval     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_evaluate_average_on_cells(const cs_xdef_t   *def,
                             cs_real_t          time_eval,
                             cs_real_t          retval[])
{
  /* Sanity checks */
  assert(def != NULL);

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    cs_evaluate_average_on_cells_by_value(def, retval);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_average_on_cells_by_analytic(def, time_eval, retval);
    break;

  case CS_XDEF_BY_ARRAY:
    cs_evaluate_average_on_cells_by_array(def, retval);

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Case not handled yet.", __func__);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EVALUATE_H__ */
