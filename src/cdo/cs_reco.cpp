/*============================================================================
 * Functions to handle the reconstruction of fields
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
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_math.h"
#include "cs_scheme_geometry.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_reco.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */

#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (e, c) --> array which can be scanned by the c2e
 *        adjacency).
 *        Case of scalar-valued array.
 *
 * \param[in] c_id       cell id
 * \param[in] c2e        cell -> edges connectivity
 * \param[in] cdoq       pointer to the additional quantities struct.
 * \param[in] array      pointer to the array of values at (e,c)
 *
 * \return the reconstructed values at the cell center
 */
/*----------------------------------------------------------------------------*/

static inline double
_scalar_ebyc2c(cs_lnum_t                    c_id,
               const cs_adjacency_t        *c2e,
               const cs_cdo_quantities_t   *cdoq,
               const double                *array)
{
  double reco_sum = 0;
  for (cs_lnum_t je = c2e->idx[c_id]; je < c2e->idx[c_id+1]; je++)
    reco_sum += cdoq->pvol_ec[je] * array[je];

  return reco_sum/cdoq->cell_vol[c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (v, c) --> array which can be scanned by the c2v
 *        adjacency).
 *        Case of scalar-valued array.
 *
 * \param[in] c_id       cell id
 * \param[in] c2v        cell -> vertices connectivity
 * \param[in] cdoq       pointer to the additional quantities struct.
 * \param[in] array      pointer to the array of values at (v,c)
 *
 * \return the reconstructed values at the cell center
 */
/*----------------------------------------------------------------------------*/

static inline double
_scalar_vbyc2c(cs_lnum_t                    c_id,
               const cs_adjacency_t        *c2v,
               const cs_cdo_quantities_t   *cdoq,
               const double                *array)
{
  double reco_sum = 0;
  for (cs_lnum_t jv = c2v->idx[c_id]; jv < c2v->idx[c_id+1]; jv++)
    reco_sum += cdoq->pvol_vc[jv] * array[jv];

  return reco_sum/cdoq->cell_vol[c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at a cell center from an array of values
 *        defined on primal vertices.
 *        Case of scalar-valued arrays.
 *
 * \param[in] c_id       cell id
 * \param[in] c2v        cell -> vertices connectivity
 * \param[in] cdoq       pointer to the additional quantities struct.
 * \param[in] array      pointer to the array of values at vertices
 *
 * \return the reconstructed values at the cell center
 */
/*----------------------------------------------------------------------------*/

static inline double
_scalar_v2c(cs_lnum_t                    c_id,
            const cs_adjacency_t        *c2v,
            const cs_cdo_quantities_t   *cdoq,
            const double                *array)
{
  double reco_sum = 0;
  for (cs_lnum_t jv = c2v->idx[c_id]; jv < c2v->idx[c_id+1]; jv++)
    reco_sum += cdoq->pvol_vc[jv] * array[c2v->ids[jv]];

  return reco_sum/cdoq->cell_vol[c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at a face center from an array of values
 *        defined on primal vertices.
 *        Case of scalar-valued arrays.
 *
 * \param[in]      f_id       face id
 * \param[in]      f2e        face -> edges connectivity
 * \param[in]      e2v        edge -> vertices connectivity
 * \param[in]      cdoq       pointer to the additional quantities struct.
 * \param[in]      array      pointer to the array of values at vertices
 *
 * \return the reconstructed values at the face center
 */
/*----------------------------------------------------------------------------*/

static double
_scalar_v2f(cs_lnum_t                    f_id,
            const cs_adjacency_t        *f2e,
            const cs_adjacency_t        *e2v,
            const cs_cdo_quantities_t   *cdoq,
            const double                *array)
{
  const cs_real_t  *xf = cs_quant_get_face_center(f_id, cdoq);

  double reco_sum = 0;
  for (cs_lnum_t i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

    const cs_lnum_t  *v_id = e2v->ids + 2*f2e->ids[i]; /* + 2*e_id */
    const cs_real_t  tef = cs_math_surftri(cdoq->vtx_coord + 3*v_id[0],
                                           cdoq->vtx_coord + 3*v_id[1],
                                           xf);

    reco_sum += (array[v_id[0]] + array[v_id[1]]) * tef;

  } /* Loop on face edges */

  return 0.5*reco_sum/cs_quant_get_face_surf(f_id, cdoq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at a cell center from an array of values
 *        defined on primal vertices.
 *        Case of vector-valued arrays.
 *
 * \param[in]      c_id       cell id
 * \param[in]      c2v        cell -> vertices connectivity
 * \param[in]      cdoq       pointer to the additional quantities struct.
 * \param[in]      array      pointer to the array of values at vertices
 * \param[in, out] reco       reconstructed values at the cell center
 */
/*----------------------------------------------------------------------------*/

static void
_vector_v2c(cs_lnum_t                    c_id,
            const cs_adjacency_t        *c2v,
            const cs_cdo_quantities_t   *cdoq,
            const double                *array,
            cs_real_t                    reco[3])
{
  reco[0] = 0, reco[1] = 0, reco[2] = 0;
  for (cs_lnum_t jv = c2v->idx[c_id]; jv < c2v->idx[c_id+1]; jv++) {

    const double  *_array = array + 3*c2v->ids[jv];
    const cs_real_t  vc_vol = cdoq->pvol_vc[jv];

    reco[0] += vc_vol * _array[0];
    reco[1] += vc_vol * _array[1];
    reco[2] += vc_vol * _array[2];

  } /* Loop on cell vertices */

  const double  invvol = 1/cdoq->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    reco[k] *= invvol;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply 1/|dual_vol| to a synchronized array of DoF vertices
 *         Parallel synchronization is done inside this function:
 *         1) A parallel sum reduction is done on the vtx_values
 *         2) Apply 1/|dual_vol| to each entry
 *
 *  \param[in]      connect    pointer to additional connectivities for CDO
 *  \param[in]      quant      pointer to additional quantities for CDO
 *  \param[in]      stride     number of entries for each vertex
 *  \param[in]      interlace  interlaced array (useful if the stride > 1)
 *  \param[in, out] array      array of DoF values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dual_vol_weight_reduction(const cs_cdo_connect_t       *connect,
                                  const cs_cdo_quantities_t    *quant,
                                  int                           stride,
                                  bool                          interlace,
                                  cs_real_t                    *array)
{
  if (array == nullptr)
    return;

  assert(stride > 0);

  const cs_lnum_t  n_vertices = quant->n_vertices;

  /* Synchronization of values at vertices */

  if (connect->vtx_ifs != nullptr)
    cs_interface_set_sum(connect->vtx_ifs,
                         n_vertices,
                         stride,
                         interlace,
                         CS_REAL_TYPE,
                         array);

  /* Compute or retrieve the dual volumes */

  cs_real_t        *build_dual_vol = nullptr;
  const cs_real_t  *dual_vol;

  if (quant->dual_vol == nullptr) {

    cs_cdo_quantities_compute_dual_volumes(quant, connect, &build_dual_vol);
    dual_vol = build_dual_vol;
  }
  else
    dual_vol = quant->dual_vol;

  /* Apply the weighting equal to 1/|dual_vol| */

  if (stride == 1) {

#pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      array[v_id] /= dual_vol[v_id];
  }
  else {

    if (interlace) {

#pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
        const cs_real_t invvol = 1. / dual_vol[v_id];
        for (int k = 0; k < stride; k++)
          array[stride * v_id + k] *= invvol;
      }
    }
    else { /* not interlace */

#pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
        const cs_real_t invvol = 1. / dual_vol[v_id];
        for (int k = 0; k < stride; k++)
          array[k * n_vertices + v_id] *= invvol;
      }

    } /* interlace ? */

  } /* stride > 1 */

  if (build_dual_vol != nullptr)
    BFT_FREE(build_dual_vol);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]   c2f           cell -> faces connectivity
 * \param[in]   cdoq          pointer to the additional quantities struct.
 * \param[in]   i_face_vals   array of DoF values for interior faces
 * \param[in]   b_face_vals   array of DoF values for border faces
 * \param[out]  cell_reco     vector-valued reconstruction inside cells
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vectors_by_ib_face_dofs(const cs_adjacency_t      *c2f,
                                     const cs_cdo_quantities_t *cdoq,
                                     const cs_real_t            i_face_vals[],
                                     const cs_real_t            b_face_vals[],
                                     cs_real_t                 *cell_reco)
{
  assert(c2f != nullptr && i_face_vals != nullptr && b_face_vals != nullptr);

  cs_array_real_fill_zero(3 * cdoq->n_cells, cell_reco);

#pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    cs_real_t *cval = cell_reco + 3 * c_id;

    /* Loop on cell faces */

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id + 1]; j++) {

      const cs_lnum_t  f_id       = c2f->ids[j];
      const cs_lnum_t  bf_id      = f_id - cdoq->n_i_faces;
      const cs_real_t *dedge_vect = cdoq->dedge_vector + 3 * j;

      if (bf_id > -1) { /* Border face */

        for (int k = 0; k < 3; k++)
          cval[k] += b_face_vals[bf_id] * dedge_vect[k];
      }
      else { /* Interior face */

        for (int k = 0; k < 3; k++)
          cval[k] += i_face_vals[f_id] * dedge_vect[k];
      }

    } /* Loop on cell faces */

    const cs_real_t invvol = 1. / cdoq->cell_vol[c_id];
    for (int k = 0; k < 3; k++)
      cval[k] *= invvol;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside a cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]  c_id          id of the cell to handle
 * \param[in]  c2f           cell -> faces connectivity
 * \param[in]  cdoq          pointer to the additional quantities struct.
 * \param[in]  face_dofs     array of DoF values at faces
 * \param[in]  local_input   true means that face_dofs is of size n_cell_faces
 * \param[out] cell_reco     vector-valued reconstruction inside cells. This
 *                           quantity should have been allocated before calling
 *                           this function
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vector_by_face_dofs(cs_lnum_t                  c_id,
                                 const cs_adjacency_t      *c2f,
                                 const cs_cdo_quantities_t *cdoq,
                                 const cs_real_t            face_dofs[],
                                 bool                       local_input,
                                 cs_real_t                 *cell_reco)
{
  assert(c2f != nullptr && cdoq != nullptr && face_dofs != nullptr);

  /* Initialization */

  cell_reco[0] = cell_reco[1] = cell_reco[2] = 0;

  /* Loop on cell faces */

  if (local_input) {

    const cs_lnum_t  s = c2f->idx[c_id], e = c2f->idx[c_id+1];
    const cs_real_t  *_dedge_vector = cdoq->dedge_vector + 3*s;

    for (cs_lnum_t j = 0; j < e-s; j++) {

      const cs_real_t  *dedge_vect = _dedge_vector + 3*j;
      for (int k = 0; k < 3; k++)
        cell_reco[k] += face_dofs[j] * dedge_vect[k];

    } /* Loop on cell faces */

  }
  else { /* face_dofs is accessed with f_id */

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->ids[j];
      const cs_real_t  *dedge_vect = cdoq->dedge_vector + 3*j;

      for (int k = 0; k < 3; k++)
        cell_reco[k] += face_dofs[f_id] * dedge_vect[k];

    } /* Loop on cell faces */

  }

  const cs_real_t  invvol = 1./cdoq->cell_vol[c_id];
  for (int k =0; k < 3; k++) cell_reco[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]  c2f          cell -> faces connectivity
 * \param[in]  cdoq         pointer to the additional quantities struct.
 * \param[in]  face_dofs    array of DoF values at faces
 * \param[out] cell_reco    vector-valued reconstruction inside cells. This
 *                          quantity should have been allocated before calling
 *                          this function
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vectors_by_face_dofs(const cs_adjacency_t       *c2f,
                                  const cs_cdo_quantities_t  *cdoq,
                                  const cs_real_t             face_dofs[],
                                  cs_real_t                  *cell_reco)
{
  assert(c2f != nullptr && cdoq != nullptr && face_dofs != nullptr);

  cs_array_real_fill_zero(3*cdoq->n_cells, cell_reco);

# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    cs_real_t  *cval = cell_reco + 3*c_id;

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->ids[j];
      const cs_real_t  *dedge_vect = cdoq->dedge_vector + 3*j;

      for (int k = 0; k < 3; k++)
        cval[k] += face_dofs[f_id] * dedge_vect[k];

    } /* Loop on cell faces */

    const cs_real_t  invvol = 1./cdoq->cell_vol[c_id];
    for (int k =0; k < 3; k++) cval[k] *= invvol;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        at primal vertices.
 *        Case of scalar-valued array.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or nullptr
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at vertices
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_v2c(cs_lnum_t                    n_cells,
                   const cs_lnum_t             *cell_ids,
                   const cs_adjacency_t        *c2v,
                   const cs_cdo_quantities_t   *cdoq,
                   const double                *array,
                   bool                         dense_ouput,
                   cs_real_t                   *reco)
{
  if (array == nullptr)
    return;
  if (n_cells < 1)
    return;

  if (cell_ids == nullptr) { /* No indirection to apply */

    assert(n_cells == cdoq->n_cells);

#pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      reco[c_id] = _scalar_v2c(c_id, c2v, cdoq, array);
  }
  else { /* There is a list of selected cells */

    if (dense_ouput) {

#pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        reco[i] = _scalar_v2c(cell_ids[i], c2v, cdoq, array);
    }
    else {

#pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++) {

        const cs_lnum_t c_id = cell_ids[i];
        reco[c_id]           = _scalar_v2c(c_id, c2v, cdoq, array);
      }

    } /* dense_ouput ? */

  } /* cell_ids != nullptr ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (v, c) --> array which can be scanned by the c2v
 *        adjacency).
 *        Case of scalar-valued array.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or nullptr
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at (v,c)
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_vbyc2c(cs_lnum_t                  n_cells,
                      const cs_lnum_t           *cell_ids,
                      const cs_adjacency_t      *c2v,
                      const cs_cdo_quantities_t *cdoq,
                      const double              *array,
                      bool                       dense_ouput,
                      cs_real_t                 *reco)
{
  if (array == nullptr)
    return;

  if (cell_ids == nullptr) { /* No indirection to apply */

    assert(c2v != nullptr && cdoq != nullptr);
    assert(n_cells == cdoq->n_cells);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      reco[c_id] = _scalar_vbyc2c(c_id, c2v, cdoq, array);
  }
  else { /* There is a list of selected cells */

    if (dense_ouput) {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        reco[i] = _scalar_vbyc2c(cell_ids[i], c2v, cdoq, array);

    }
    else {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++) {

        const cs_lnum_t  c_id = cell_ids[i];
        reco[c_id] = _scalar_vbyc2c(c_id, c2v, cdoq, array);

      }

    } /* dense_ouput ? */

  } /* cell_ids != nullptr ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (e, c) --> array which can be scanned by the c2e
 *        adjacency).
 *        Case of scalar-valued array.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or nullptr
 * \param[in]      c2e           cell -> edges connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at (e,c)
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_ebyc2c(cs_lnum_t                    n_cells,
                      const cs_lnum_t             *cell_ids,
                      const cs_adjacency_t        *c2e,
                      const cs_cdo_quantities_t   *cdoq,
                      const double                *array,
                      bool                         dense_ouput,
                      cs_real_t                   *reco)
{
  if (array == nullptr)
    return;

  assert(cdoq->pvol_ec != nullptr);

  if (cell_ids == nullptr) { /* No indirection to apply */

    assert(c2e != nullptr && cdoq != nullptr);
    assert(n_cells == cdoq->n_cells);

#pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      reco[c_id] = _scalar_ebyc2c(c_id, c2e, cdoq, array);
  }
  else { /* There is a list of selected cells */

    if (dense_ouput) {

#pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        reco[i] = _scalar_ebyc2c(cell_ids[i], c2e, cdoq, array);
    }
    else {

#pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++) {

        const cs_lnum_t c_id = cell_ids[i];
        reco[c_id]           = _scalar_ebyc2c(c_id, c2e, cdoq, array);
      }

    } /* dense_ouput ? */

  } /* cell_ids != nullptr ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at all cell centers from an array of values
 *        defined on primal vertices.
 *        Case of vector-valued fields.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or nullptr
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at vertices
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vector_v2c(cs_lnum_t                  n_cells,
                   const cs_lnum_t           *cell_ids,
                   const cs_adjacency_t      *c2v,
                   const cs_cdo_quantities_t *cdoq,
                   const double              *array,
                   bool                       dense_ouput,
                   cs_real_t                 *reco)
{
  if (array == nullptr)
    return;

  assert(c2v != nullptr && cdoq != nullptr);

  if (cell_ids == nullptr) { /* No indirection to apply */

    assert(c2v != nullptr && cdoq != nullptr);
    assert(n_cells == cdoq->n_cells);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      _vector_v2c(c_id, c2v, cdoq, array, reco + 3 * c_id);
  }
  else { /* There is a list of selected cells */

    if (dense_ouput) {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        _vector_v2c(cell_ids[i], c2v, cdoq, array, reco + 3*i);

    }
    else {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++) {

        const cs_lnum_t  c_id = cell_ids[i];
        _vector_v2c(c_id, c2v, cdoq, array, reco + 3*c_id);

      } /* Loop on selected cells */

    } /* dense_ouput ? */

  } /* cell_ids != nullptr ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at the face center from an array of values
 *        defined on primal vertices.
 *        Case of scalar-valued arrays.
 *
 * \param[in]      n_faces      number of faces
 * \param[in]      face_ids     list of face ids (interior/border faces) or
 * nullptr \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to the additional quantities struct.
 * \param[in]      array        pointer to the array of values at vertices
 * \param[in]      dense_ouput  apply cell_ids on the reco array
 * \param[in, out] reco         value of the reconstruction at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_v2f(cs_lnum_t                     n_faces,
                   const cs_lnum_t              *face_ids,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *cdoq,
                   const double                 *array,
                   bool                          dense_ouput,
                   cs_real_t                    *reco)
{
  if (array == nullptr)
    return;

  const cs_adjacency_t  *f2e = connect->f2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  if (face_ids == nullptr) { /* No indirection to apply */

#pragma omp parallel for if (n_faces > CS_THR_MIN)
    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
      reco[f_id] = _scalar_v2f(f_id, f2e, e2v, cdoq, array);
  }
  else {

    if (dense_ouput) {

#pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_faces; i++)
        reco[i] = _scalar_v2f(face_ids[i], f2e, e2v, cdoq, array);
    }
    else {

#pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_faces; i++) {

        const cs_lnum_t f_id = face_ids[i];
        reco[f_id]           = _scalar_v2f(f_id, f2e, e2v, cdoq, array);

      } /* Loop on selected faces */

    } /* dense_ouput ? */

  } /* face_ids != nullptr ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct at cell centers and face centers a vertex-based field
 *        Linear interpolation. If p_reco_c and/or p_reco_f are not allocated,
 *        this is done in this subroutine.
 *        Case of scalar-valued arrays.
 *
 *  \param[in]      connect   pointer to additional connectivities for CDO
 *  \param[in]      cdoq      pointer to additional quantities for CDO
 *  \param[in]      dof       pointer to the field of vtx-based DoFs
 *  \param[in, out] p_reco_c  reconstructed values at cell centers
 *  \param[in, out] p_reco_f  reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_v2c_v2f(const cs_cdo_connect_t    *connect,
                       const cs_cdo_quantities_t *cdoq,
                       const double              *dof,
                       double                    *p_reco_c[],
                       double                    *p_reco_f[])
{
  if (dof == nullptr)
    return;

  double *crec = *p_reco_c, *frec = *p_reco_f;

  /* Allocate arrays if necessary */

  if (crec == nullptr)
    BFT_MALLOC(crec, cdoq->n_cells, double);
  if (frec == nullptr)
    BFT_MALLOC(frec, cdoq->n_faces, double);

  const cs_adjacency_t *c2v = connect->c2v;
  const cs_adjacency_t *f2e = connect->f2e;
  const cs_adjacency_t *e2v = connect->e2v;

  /* Reconstruction at cell centers */

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++)
    crec[c_id] = _scalar_v2c(c_id, c2v, cdoq, dof);

  /* Reconstruction at face centers */

  for (cs_lnum_t f_id = 0; f_id < cdoq->n_faces; f_id++)
    frec[f_id] = _scalar_v2f(f_id, f2e, e2v, cdoq, dof);

  /* Return pointers */

  *p_reco_c = crec;
  *p_reco_f = frec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct at face centers by a cell-based field
 *        weighted average.
 *
 *  \param[in]      connect   pointer to additional connectivities for CDO
 *  \param[in]      cdoq      pointer to additional quantities for CDO
 *  \param[in]      p_c       dofs at cell centers
 *  \param[in, out] p_reco_f  reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_c2f(const cs_cdo_connect_t     *connect,
                   const cs_cdo_quantities_t  *cdoq,
                   const double               *p_c,
                   double                     *p_reco_f)
{
  if (p_c == NULL || p_reco_f == NULL)
    return;

  /* Allocate arrays if necessary */

  const cs_adjacency_t  *c2f = connect->c2f;

  assert(cdoq->pvol_fc != NULL);

  const double *pvol_fc = cdoq->pvol_fc;

  double *pvol_f;
  BFT_MALLOC(pvol_f, cdoq->n_faces, double);

  memset(pvol_f, 0.0, cdoq->n_faces*sizeof(double));
  memset(p_reco_f, 0.0, cdoq->n_faces*sizeof(double));

  /* Reconstruction at face centers */

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {
      cs_lnum_t f_id = c2f->ids[j];
      pvol_f[f_id] += pvol_fc[j];
      p_reco_f[f_id] += pvol_fc[j]*p_c[c_id];
    }
  }

# pragma omp parallel for if (cdoq->n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < cdoq->n_faces; f_id++)
    p_reco_f[f_id] /= pvol_f[f_id];

  BFT_FREE(pvol_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a scalar-valued array at vertices from a scalar-valued
 *         array at cells.
 *
 * \param[in]      connect   pointer to additional connectivities for CDO
 * \param[in]      quant     pointer to the additional quantities for CDO
 * \param[in]      cell_val  array of scalar-valued values at cells
 * \param[in, out] vtx_val   array of scalar-valued values at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scal_pv_from_pc(const cs_cdo_connect_t    *connect,
                        const cs_cdo_quantities_t *quant,
                        const cs_real_t           *cell_val,
                        cs_real_t                 *vtx_val)
{
  if (cell_val == nullptr || vtx_val == nullptr)
    return;
  assert(quant != nullptr && connect != nullptr);

  cs_array_real_fill_zero(quant->n_vertices, vtx_val);

  const cs_adjacency_t *c2v = connect->c2v;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t cval = cell_val[c_id];
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id + 1]; j++)
      vtx_val[c2v->ids[j]] += quant->pvol_vc[j] * cval;

  } /* Loop on cells */

  cs_reco_dual_vol_weight_reduction(connect, quant, 1, true, vtx_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a vector-valued array at vertices from a vector-valued
 *         array at cells.
 *
 * \param[in]      connect   pointer to additional connectivities for CDO
 * \param[in]      quant     pointer to the additional quantities for CDO
 * \param[in]      cell_val  array of vector-valued values at cells
 * \param[in, out] vtx_val   array of vector-valued values at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vect_pv_from_pc(const cs_cdo_connect_t    *connect,
                        const cs_cdo_quantities_t *quant,
                        const cs_real_t           *cell_val,
                        cs_real_t                 *vtx_val)
{
  if (cell_val == nullptr || vtx_val == nullptr)
    return;
  assert(quant != nullptr && connect != nullptr);

  cs_array_real_fill_zero(3 * quant->n_vertices, vtx_val);

  const cs_adjacency_t *c2v = connect->c2v;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t *cval = cell_val + 3 * c_id;
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id + 1]; j++) {

      const cs_real_t vc_vol = quant->pvol_vc[j];
      cs_real_t      *vval   = vtx_val + 3 * c2v->ids[j];

      vval[0] += vc_vol * cval[0];
      vval[1] += vc_vol * cval[1];
      vval[2] += vc_vol * cval[2];

    } /* Loop on cell vertices */
  } /* Loop on cells */

  cs_reco_dual_vol_weight_reduction(connect, quant, 3, true, vtx_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector at the cell center from an array of
 *         values defined on dual faces lying inside each cell.
 *         This array is scanned thanks to the c2e connectivity.
 *
 *  \param[in]      c_id     cell id
 *  \param[in]      c2e      cell -> edges connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_at_cell_center(cs_lnum_t                  c_id,
                             const cs_adjacency_t      *c2e,
                             const cs_cdo_quantities_t *quant,
                             const double              *array,
                             cs_real_3_t                val_xc)
{
  /* Initialization */

  val_xc[0] = val_xc[1] = val_xc[2] = 0.;

  if (array == nullptr)
    return;

  for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id + 1]; j++) {

    const cs_real_t *e_vect = quant->edge_vector + 3 * c2e->ids[j];
    for (int k = 0; k < 3; k++)
      val_xc[k] += array[j] * e_vect[k];

  } /* Loop on cell edges */

  const double invvol = 1 / quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    val_xc[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector inside the cell c.
 *         array is scanned thanks to the c2e connectivity. Pointer is already
 *         located at the beginning of the cell sequence.
 *         Reconstruction used is based on DGA (stabilization = 1/d where d is
 *         the space dimension)
 *
 *  \param[in]      cm        pointer to a cs_cell_mesh_t structure
 *  \param[in]      array     local pointer to the array of values
 *  \param[in, out] val_c     value of the reconstructed vector in the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_in_cell(const cs_cell_mesh_t *cm,
                      const cs_real_t      *array,
                      cs_real_3_t           val_c)
{
  /* Initialization */

  val_c[0] = val_c[1] = val_c[2] = 0.;

  if (array == nullptr)
    return;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PEQ));

  const double invvol = 1 / cm->vol_c;

  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_quant_t peq          = cm->edge[e];
    const cs_real_t  edge_contrib = array[e] * peq.meas;

    for (int k = 0; k < 3; k++)
      val_c[k] += edge_contrib * peq.unitv[k];

  } /* Loop on cell edges */

  for (int k = 0; k < 3; k++)
    val_c[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector inside pec which is a volume
 *         surrounding the edge e inside the cell c.
 *         array is scanned thanks to the c2e connectivity. Pointer is already
 *         located at the beginning of the cell sequence.
 *         Reconstruction used is based on DGA (stabilization = 1/d where d is
 *         the space dimension)
 *
 *  \param[in]      cm        pointer to a cs_cell_mesh_t structure
 *  \param[in]      e         local edge id
 *  \param[in]      array     local pointer to the array of values
 *  \param[in, out] val_pec   value of the reconstruction in pec
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_in_pec(const cs_cell_mesh_t *cm,
                     short int             e,
                     const cs_real_t      *array,
                     cs_real_3_t           val_pec)
{
  /* Initialize values */

  val_pec[0] = val_pec[1] = val_pec[2] = 0.;

  if (array == nullptr)
    return;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ));

  cs_real_3_t val_c = { 0., 0., 0. };

  /* Compute val_c */

  for (short int _e = 0; _e < cm->n_ec; _e++) {

    const cs_quant_t _peq = cm->edge[_e];

    for (int k = 0; k < 3; k++)
      val_c[k] += array[_e] * _peq.meas * _peq.unitv[k];

  } /* Loop on cell edges */

  const double invvol = 1 / cm->vol_c;

  /* Compute the consistency part related to this cell */

  for (int k = 0; k < 3; k++)
    val_c[k] *= invvol;

  /* Compute the reconstruction inside pec */

  const cs_quant_t peq   = cm->edge[e];
  const cs_nvec3_t dfq   = cm->dface[e];
  const double     ecoef = (array[e] - dfq.meas * _dp3(dfq.unitv, val_c))
                       / (dfq.meas * _dp3(dfq.unitv, peq.unitv));

  for (int k = 0; k < 3; k++)
    val_pec[k] = val_c[k] + ecoef * peq.unitv[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at the cell center a field of edge-based DoFs
 *
 *  \param[in]      c_id    cell id
 *  \param[in]      c2e     cell -> edges connectivity
 *  \param[in]      quant   pointer to the additional quantities struct.
 *  \param[in]      dof     pointer to the field of edge-based DoFs
 *  \param[in, out] reco    value of the reconstructed field at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dof(cs_lnum_t                  c_id,
                      const cs_adjacency_t      *c2e,
                      const cs_cdo_quantities_t *quant,
                      const double              *dof,
                      double                     reco[])
{
  if (dof == nullptr)
    return;

  /* Initialize value */

  reco[0] = reco[1] = reco[2] = 0.0;

  const cs_lnum_t *c2e_idx = c2e->idx + c_id;
  const cs_lnum_t *c2e_ids = c2e->ids + c2e_idx[0];
  const cs_real_t *dface   = quant->dface_normal + 3 * c2e_idx[0];
  for (cs_lnum_t i = 0; i < c2e_idx[1] - c2e_idx[0]; i++) {

    /* Dual face quantities */

    const double val = dof[c2e_ids[i]]; /* Edge value */
    for (int k = 0; k < 3; k++)
      reco[k] += val * dface[3 * i + k];

  } /* End of loop on cell edges */

  /* Divide by the cell volume */

  const double invvol = 1 / quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    reco[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at each cell center a field of edge-based DoFs
 *
 *  \param[in]      connect   pointer to the connectivity struct.
 *  \param[in]      quant     pointer to the additional quantities struct.
 *  \param[in]      dof       pointer to the field of edge-based DoFs
 *  \param[in, out] p_ccrec   pointer to the reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dofs(const cs_cdo_connect_t    *connect,
                       const cs_cdo_quantities_t *quant,
                       const double              *dof,
                       double                    *p_ccrec[])
{
  assert(connect->c2e != nullptr); /* Sanity check */

  double *ccrec = *p_ccrec;

  if (dof == nullptr)
    return;

  /* Allocate reconstructed vector field at each cell barycenter */

  if (ccrec == nullptr)
    BFT_MALLOC(ccrec, 3 * quant->n_cells, double);

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (int c_id = 0; c_id < quant->n_cells; c_id++)
    cs_reco_ccen_edge_dof(c_id, connect->c2e, quant, dof, ccrec + 3 * c_id);

  /* Return pointer */

  *p_ccrec = ccrec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a cell-wise constant curl from the knowledge of the
 *         circulation at primal edges
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      circ     pointer to the array of circulations at edges
 * \param[in, out] p_curl   pointer to value of the reconstructed curl inside
 *                          cells (allocated if set to nullptr)
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_curl_by_edge_dofs(const cs_cdo_connect_t    *connect,
                               const cs_cdo_quantities_t *quant,
                               const cs_real_t           *circ,
                               cs_real_t                **p_curl)
{
  if (circ == nullptr)
    return;

  assert(connect != nullptr && quant != nullptr);

  cs_real_t *curl_vectors = *p_curl;
  if (curl_vectors == nullptr)
    BFT_MALLOC(curl_vectors, 3 * quant->n_cells, cs_real_t);

  cs_real_t *face_curl = nullptr;
  cs_cdo_connect_discrete_curl(connect, circ, &face_curl);

  cs_reco_cell_vectors_by_face_dofs(connect->c2f, quant, face_curl,
                                    curl_vectors);

  BFT_FREE(face_curl);

  /* Returns pointer */

  *p_curl = curl_vectors;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the mean-value of the gradient field with DoFs arising
 *         from a face-based scheme (values at face center and cell center)
 *         The reconstruction only deals with the consistent part so that there
 *         is no distinction between Fb schemes
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      p_c      pointer to the array of values in cells
 * \param[in]      p_f      pointer to the array of values on faces
 * \param[in, out] grd_c    value of the reconstructed gradient at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_cell_from_fb_dofs(cs_lnum_t                    c_id,
                               const cs_cdo_connect_t      *connect,
                               const cs_cdo_quantities_t   *quant,
                               const cs_real_t             *p_c,
                               const cs_real_t             *p_f,
                               cs_real_t                    grd_c[])
{
  /* Initialize the value to return */

  grd_c[0] = grd_c[1] = grd_c[2] = 0.;

  if (p_c == nullptr || p_f == nullptr)
    return;
  assert(connect != nullptr && quant != nullptr); /* sanity checks */

  const cs_real_t  _p_c = p_c[c_id];
  const cs_lnum_t  s = connect->c2f->idx[c_id], e = connect->c2f->idx[c_id+1];
  const cs_lnum_t  *c2f_ids = connect->c2f->ids + s;
  const short int  *c2f_sgn = connect->c2f->sgn + s;

  for (cs_lnum_t  i = 0; i < e-s; i++) {

    const cs_lnum_t  f_id = c2f_ids[i];
    const cs_real_t  *fq = cs_quant_get_face_vector_area(f_id, quant);
    for (int k = 0; k < 3; k++)
      grd_c[k] += c2f_sgn[i] * (p_f[f_id] - _p_c) * fq[k];

  } /* Loop on cell faces */

  const cs_real_t  invvol = 1./quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    grd_c[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the mean-value of the tensor gradient field with DoFs
 *         arising from a face-based scheme (vector-valued at face center and
 *         cell center) The reconstruction only deals with the consistent part
 *         so that there is no distinction between Fb schemes
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      u_c      pointer to the array of values in cells
 * \param[in]      u_f      pointer to the array of values on faces
 * \param[in, out] grd_c    value of the reconstructed gradient at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_33_cell_from_fb_dofs(cs_lnum_t                    c_id,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  const cs_real_t             *u_c,
                                  const cs_real_t             *u_f,
                                  cs_real_t                    grd_c[])
{
  /* Initialize the value to return */

  grd_c[0] = grd_c[1] = grd_c[2] = 0.;
  grd_c[3] = grd_c[4] = grd_c[5] = 0.;
  grd_c[6] = grd_c[7] = grd_c[8] = 0.;

  if (u_c == nullptr || u_f == nullptr)
    return;
  assert(connect != nullptr && quant != nullptr); /* sanity checks */

  const cs_real_t  *_u_c = u_c + 3*c_id;
  const cs_lnum_t  s = connect->c2f->idx[c_id], e = connect->c2f->idx[c_id+1];
  const cs_lnum_t  *c2f_ids = connect->c2f->ids + s;
  const short int  *c2f_sgn = connect->c2f->sgn + s;

  for (cs_lnum_t  i = 0; i < e-s; i++) {

    const cs_lnum_t  f_id = c2f_ids[i];
    const cs_real_t  *fq = cs_quant_get_face_vector_area(f_id, quant);
    const cs_real_t  *_u_f = u_f + 3*f_id;

    for (int ki = 0; ki < 3; ki++) {
      const cs_real_t  gki = c2f_sgn[i] * (_u_f[ki] - _u_c[ki]);
      for (int kj = 0; kj < 3; kj++)
        grd_c[3*ki+kj] += gki * fq[kj];
    }

  } /* Loop on cell faces */

  const cs_real_t  invvol = 1./quant->cell_vol[c_id];
  for (int k = 0; k < 9; k++)
    grd_c[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the constant gradient vector in a cell (the mean value)
 *        from the value at mesh vertices.
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      pdi      pointer to the array of values
 * \param[in, out] grdc     value of the reconstructed gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_cell_from_pv(cs_lnum_t                    c_id,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_real_t             *pdi,
                          cs_real_t                    grdc[])
{
  grdc[0] = grdc[1] = grdc[2] = 0.;

  if (pdi == nullptr)
    return;

  const cs_adjacency_t  *e2v = connect->e2v;
  const cs_adjacency_t  *c2e = connect->c2e;
  const cs_lnum_t  *c2e_idx = c2e->idx + c_id;
  const cs_lnum_t  *c2e_ids = c2e->ids + c2e_idx[0];
  const cs_real_t  *dface = quant->dface_normal + 3*c2e_idx[0];

  for (cs_lnum_t i = 0; i < c2e_idx[1] - c2e_idx[0]; i++) {

    const cs_lnum_t  shift_e = 2*c2e_ids[i];
    const short int  sgn_va = e2v->sgn[shift_e];
    const cs_lnum_t  va = e2v->ids[shift_e], vb = e2v->ids[shift_e+1];

    const cs_real_t  gdi_e = sgn_va*(pdi[va] - pdi[vb]);

    for (int k = 0; k < 3; k++)
      grdc[k] += gdi_e * dface[3*i+k];

  } /* Loop on cell edges */

  /* Divide by the cell volume */

  const double  invvol = 1/quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    grdc[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the constant gradient vector in a cell (the mean value)
 *        from the value at mesh vertices. Case of two scalar fields.
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      p1di     pointer to the array of values
 * \param[in]      p2di     pointer to the array of values
 * \param[in, out] grd1c    value of the reconstructed gradient for p1
 * \param[in, out] grd2c    value of the reconstructed gradient for p2
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_2grad_cell_from_pv(cs_lnum_t                    c_id,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_real_t             *p1di,
                           const cs_real_t             *p2di,
                           cs_real_t                    grd1c[],
                           cs_real_t                    grd2c[])
{
  grd1c[0] = grd1c[1] = grd1c[2] = 0.;
  grd2c[0] = grd2c[1] = grd2c[2] = 0.;

  if (p1di == nullptr || p2di == nullptr)
    return;

  const cs_adjacency_t  *e2v = connect->e2v;
  const cs_adjacency_t  *c2e = connect->c2e;
  const cs_lnum_t  *c2e_idx = c2e->idx + c_id;
  const cs_lnum_t  *c2e_ids = c2e->ids + c2e_idx[0];
  const cs_real_t  *dface = quant->dface_normal + 3*c2e_idx[0];

  for (cs_lnum_t i = 0; i < c2e_idx[1] - c2e_idx[0]; i++) {

    const cs_lnum_t  shift_e = 2*c2e_ids[i];
    const short int  sgn_va = e2v->sgn[shift_e];
    const cs_lnum_t  va = e2v->ids[shift_e], vb = e2v->ids[shift_e+1];

    const cs_real_t  g1di_e = sgn_va*(p1di[va] - p1di[vb]);
    const cs_real_t  g2di_e = sgn_va*(p2di[va] - p2di[vb]);

    for (int k = 0; k < 3; k++) {
      grd1c[k] += g1di_e * dface[3*i+k];
      grd2c[k] += g2di_e * dface[3*i+k];
    }

  } /* Loop on cell edges */

  /* Divide by the cell volume */

  const double  invvol = 1/quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++) {
    grd1c[k] *= invvol;
    grd2c[k] *= invvol;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the
 *        given flux array. This array is stored in the same order as cm->f_ids
 *        Scalar-valued face DoFs are related to the normal flux across primal
 *        faces.
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes         array of normal fluxes on primal faces
 * \param[out]     cell_reco      vector-valued reconstruction inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_vect_from_flux(const cs_cell_mesh_t    *cm,
                               const cs_real_t         *fluxes,
                               cs_real_t               *cell_reco)
{
  if (fluxes == nullptr)
    return;

  assert(cell_reco != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialization */

  cell_reco[0] = cell_reco[1] = cell_reco[2] = 0.;

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    /* Related dual edge quantity */

    const cs_nvec3_t  deq = cm->dedge[f];
    const cs_real_t  coef = fluxes[f] * deq.meas;

    cell_reco[0] += coef*deq.unitv[0];
    cell_reco[1] += coef*deq.unitv[1];
    cell_reco[2] += coef*deq.unitv[2];

  } /* Loop on cell faces */

  const cs_real_t  invvol = 1./cm->vol_c;
  for (int k =0; k < 3; k++) cell_reco[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      i_face_vals    array of DoF values for interior faces
 * \param[in]      b_face_vals    array of DoF values for border faces
 * \param[out]     cell_reco      vector-valued reconstruction inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_vect_from_face_dofs(const cs_cell_mesh_t    *cm,
                                    const cs_real_t          i_face_vals[],
                                    const cs_real_t          b_face_vals[],
                                    cs_real_t               *cell_reco)
{
  assert(cm != nullptr && i_face_vals != nullptr && b_face_vals != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialization */

  cell_reco[0] = cell_reco[1] = cell_reco[2] = 0.;

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    /* Related dual edge quantity */

    const cs_lnum_t  f_id = cm->f_ids[f];
    const cs_nvec3_t  deq = cm->dedge[f];
    const cs_real_t  coef = (cm->f_ids[f] < cm->bface_shift) ?
      i_face_vals[f_id] * deq.meas :
      b_face_vals[f_id - cm->bface_shift] * deq.meas;

    cell_reco[0] += coef*deq.unitv[0];
    cell_reco[1] += coef*deq.unitv[1];
    cell_reco[2] += coef*deq.unitv[2];

  } /* Loop on cell faces */

  const cs_real_t  invvol = 1./cm->vol_c;
  for (int k =0; k < 3; k++) cell_reco[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at a point inside a cell
 *         The scalar potential has DoFs located at primal vertices
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      pdi            array of DoF values at vertices
 * \param[out]     cell_gradient  gradient inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_grad_from_scalar_pv(const cs_cell_mesh_t    *cm,
                                    const cs_real_t          pdi[],
                                    cs_real_t               *cell_gradient)
{
  assert(cm != nullptr && pdi != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PVQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_DFQ));

  /* Reconstruct a constant gradient inside the current cell */

  cell_gradient[0] = cell_gradient[1] = cell_gradient[2] = 0;
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_lnum_t  v0 = cm->v_ids[cm->e2v_ids[2*e]];
    const cs_lnum_t  v1 = cm->v_ids[cm->e2v_ids[2*e+1]];
    const cs_real_t  coeff = cm->e2v_sgn[e]*(pdi[v0]-pdi[v1])*cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      cell_gradient[k] += coeff * cm->dface[e].unitv[k];

  } /* Loop on cell edges */

  const cs_real_t  invcell = 1./cm->vol_c;
  for (int k = 0; k < 3; k++) cell_gradient[k] *= invcell;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at a point inside a cell
 *         The scalar potential has DoFs located at primal vertices
 *
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      pdi          array of DoF values at vertices
 * \param[in]      length_xcxp  lenght of the segment [x_c, x_p]
 * \param[in]      unitv_xcxp   unitary vector pointed from x_c to x_p
 * \param[in, out] wbuf         pointer to a temporary buffer
 *
 * \return the value of the reconstructed potential at the evaluation point
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_reco_cw_scalar_pv_inside_cell(const cs_cell_mesh_t    *cm,
                                 const cs_real_t          pdi[],
                                 const cs_real_t          length_xcxp,
                                 const cs_real_t          unitv_xcxp[],
                                 cs_real_t                wbuf[])
{
  assert(cm != nullptr && pdi != nullptr && wbuf != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PVQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_DFQ));

  cs_real_t  *_pv = wbuf;  /* Local value of the potential field */

  /* Reconstruct the value at the cell center */

  cs_real_t  pc = 0.;
  for (short int v = 0; v < cm->n_vc; v++) {
    _pv[v] = pdi[cm->v_ids[v]];
    pc += cm->wvc[v] * _pv[v];
  }

  /* Reconstruct a constant gradient inside the current cell */

  cs_real_3_t  gcell = {0., 0., 0.};
  for (short int e = 0; e < cm->n_ec; e++) {
    const cs_real_t  ge =
      cm->e2v_sgn[e]*( _pv[cm->e2v_ids[2*e]] - _pv[cm->e2v_ids[2*e+1]]);
    const cs_real_t  coef_e = ge * cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      gcell[k] += coef_e * cm->dface[e].unitv[k];
  }

  const cs_real_t  invcell = 1./cm->vol_c;
  for (int k = 0; k < 3; k++) gcell[k] *= invcell;

  /* Evaluation at the given point */

  cs_real_t  p_rec = pc;
  p_rec += length_xcxp * cs_math_3_dot_product(gcell, unitv_xcxp);

  return  p_rec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the weighted (by volume) gradient inside a given primal
 *        cell for the related vertices.
 *        Use the WBS algo. for approximating the gradient.
 *        The computation takes into account a subdivision into tetrahedra of
 *        the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices + cell
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] vgrd     gradient at vertices inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_vgrd_wbs_from_pvc(const cs_cell_mesh_t   *cm,
                             const cs_real_t        *pot,
                             cs_cell_builder_t      *cb,
                             cs_real_t              *vgrd)
{
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  | CS_FLAG_COMP_HFQ));

  constexpr cs_real_t c_1ov3 = 1./3.;

  cs_real_3_t  grd_c, grd_v1, grd_v2;

  /* Temporary buffers */

  cs_real_3_t  *u_vc = cb->vectors;
  double  *l_vc = cb->values;

  const double  *p_v = pot;
  const double  p_c = pot[cm->n_vc];

  /* Reset local fluxes */

  for (int i = 0; i < 3*cm->n_vc; i++) vgrd[i] = 0.;

  /* Store segments xv --> xc for this cell */

  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute geometrical quantities for the current face:
       - the weighting of each triangle defined by a base e and an apex f
       - volume of each sub-tetrahedron pef_c
       - the gradient of the Lagrange function related xc in p_{f,c} */

    const cs_real_t  ohf = -cm->f_sgn[f]/cm->hfc[f];
    for (int k = 0; k < 3; k++) grd_c[k] = ohf * pfq.unitv[k];

    /* Compute the reconstructed value of the potential at p_f */

    double  p_f = 0.;
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];

      p_f += cm->tef[i]*(  p_v[e2v[0]] + p_v[e2v[1]] );

    }
    p_f *= 0.5/pfq.meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */

    const cs_real_t  hf_coef = c_1ov3 * cm->hfc[f];
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
      const short int  v1 = e2v[0], v2 = e2v[1];

      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                        grd_v1, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c} */

      const cs_real_t  sefv_vol = 0.5 * hf_coef * cm->tef[i]; /* 0.5 pef_vol */
      cs_real_t  _grd;
      for (int k = 0; k < 3; k++) {
        _grd = sefv_vol * ( dp_cf           * grd_c[k]  +
                            (p_v[v1] - p_f) * grd_v1[k] +
                            (p_v[v2] - p_f) * grd_v2[k]);
        vgrd[3*v1 + k] += _grd;
        vgrd[3*v2 + k] += _grd;
      }

    } /* Loop on face edges */

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the mean value of a gradient inside a given primal cell.
 *          Use the WBS algo. for approximating the gradient.
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices + cell
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] cgrd     mean-value of the cell gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cgrd_wbs_from_pvc(const cs_cell_mesh_t   *cm,
                             const cs_real_t        *pot,
                             cs_cell_builder_t      *cb,
                             cs_real_t              *cgrd)
{
  constexpr cs_real_t c_1ov3 = 1./3.;

  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  | CS_FLAG_COMP_HFQ));

  cs_real_3_t  grd_c, grd_v1, grd_v2;

  /* Temporary buffers */

  cs_real_3_t  *u_vc = cb->vectors;
  double  *l_vc = cb->values;

  cgrd[0] = cgrd[1] = cgrd[2] = 0.;

  const double  *p_v = pot;
  const double  p_c = pot[cm->n_vc];

  /* Store segments xv --> xc for this cell */

  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute geometrical quantities for the current face:
       - the weighting of each triangle defined by a base e and an apex f
       - volume of each sub-tetrahedron pef_c
       - the gradient of the Lagrange function related xc in p_{f,c} */

    const cs_real_t  ohf = -cm->f_sgn[f]/cm->hfc[f];
    for (int k = 0; k < 3; k++) grd_c[k] = ohf * pfq.unitv[k];

    /* Compute the reconstructed value of the potential at p_f */

    double  p_f = 0.;
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];

      p_f += cm->tef[i]*( p_v[e2v[0]] + p_v[e2v[1]] );

    }
    p_f *= 0.5/pfq.meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */

    const cs_real_t  hf_coef = c_1ov3 * cm->hfc[f];
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
      const short int  v1 = e2v[0], v2 = e2v[1];

      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                        grd_v1, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c} */

      const cs_real_t  pefv_vol = hf_coef * cm->tef[i];
      for (int k = 0; k < 3; k++)
        cgrd[k] += pefv_vol * ( dp_cf          * grd_c[k]  +
                               (p_v[v1] - p_f) * grd_v1[k] +
                               (p_v[v2] - p_f) * grd_v2[k]);

    } /* Loop on face edges */

  } /* Loop on cell faces */

  const double  invvol = 1/cm->vol_c;
  for (int k = 0; k < 3; k++) cgrd[k] *= invvol;
}

/*----------------------------------------------------------------------------*/

#undef _dp3
END_C_DECLS
