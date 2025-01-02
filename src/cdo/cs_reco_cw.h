#ifndef __CS_RECO_CW_H__
#define __CS_RECO_CW_H__

/*============================================================================
 * Functions to handle the cell-wise reconstruction of fields relying on the
 * cell mesh structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cdo/cs_cdo_local.h"

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
 * \brief Reconstruct a scalar value at the cell center.
 *        array is assumed to be interlaced and of size stride*n_vertices
 *
 * \param[in]      stride    number of values for each vertex
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      array     array of values
 * \param[in, out] reco      reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_stride_v2c(int                        stride,
                      const cs_cell_mesh_t      *cm,
                      const cs_real_t           *array,
                      cs_real_t                 *reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct a scalar value at the cell center.
 *        array is scanned thanks to the c2v connectivity. Pointer is already
 *        located at the beginning of the cell sequence, i.e. a shift equal to
 *        stride*c2v->idx[cm->c_id] has been done.
 *        array is assumed to be interlaced
 *
 * \param[in]      stride    number of values for each vertex
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      array     array of values for each couple (v,c)
 * \param[in, out] reco      reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_stride_vbyc2c(int                        stride,
                         const cs_cell_mesh_t      *cm,
                         const cs_real_t           *array,
                         cs_real_t                 *reco);

/*============================================================================
 * Inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at the face center from an array of values
 *        defined on primal vertices attached to a face.
 *
 * \param[in] fm     pointer to cs_face_mesh_t structure
 * \param[in] p_v    pointer to an array of values (local to this face)
 *
 * \return the value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_v2f_fw(const cs_face_mesh_t    *fm,
                         const cs_real_t         *p_v)
{
  cs_real_t  p_f = 0.;

  if (p_v == NULL)
    return p_f;

  const cs_quant_t  pfq = fm->face;

  for (short int e = 0; e < fm->n_ef; e++)
    p_f += (p_v[fm->e2v_ids[2*e]] + p_v[fm->e2v_ids[2*e+1]]) * fm->tef[e];
  p_f *= 0.5 / pfq.meas;

  return p_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at the face center from an array of values
 *        defined on primal vertices.
 *
 * \param[in] f      id of the face in the cellwise numbering
 * \param[in] cm     pointer to cs_cell_mesh_t structure
 * \param[in] p_v    pointer to an array of values (local to the cell)
 *
 * \return the value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_v2f_loc(const short int          f,
                          const cs_cell_mesh_t    *cm,
                          const cs_real_t         *p_v)
{
  cs_real_t  p_f = 0.;

  if (p_v == NULL)
    return p_f;

  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ |
                       CS_FLAG_COMP_FE));

  for (int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++) {
    const short int  *v = cm->e2v_ids + 2*cm->f2e_ids[ie];
    p_f += (p_v[v[0]] + p_v[v[1]]) * cm->tef[ie];
  }
  p_f *= 0.5 / cm->face[f].meas;

  return p_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value of a scalar potential at the cell center from
 *        an array of values defined on primal vertices.
 *        Algorithm based on the cs_cell_mesh_t structure.
 *
 * \param[in] cm    pointer to a cs_cell_mesh_t structure
 * \param[in] p_v   pointer to the array of values at vertices (size: n_vc)
 *
 * \return the value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_v2c_loc(const cs_cell_mesh_t     *cm,
                          const cs_real_t          *p_v)
{
  cs_real_t  p_c = 0.;

  if (p_v == NULL || cm == NULL)
    return p_c;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ)); /* Sanity check */

  /* Reconstruct the value at the cell center */

  for (short int v = 0; v < cm->n_vc; v++)
    p_c += cm->wvc[v] * p_v[v];

  return p_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct a scalar value at the cell center from an array of values
 *        located at vertices (for the full mesh)
 *
 * \param[in] cm         pointer to a cs_cell_mesh_t structure
 * \param[in] array      array of values at vertices (size:n_vertices)
 *
 * \return the reconstructed value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_v2c(const cs_cell_mesh_t      *cm,
                      const cs_real_t           *array)
{
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  cs_real_t  val_c = 0;

  if (array == NULL)
    return val_c;

  for (short int v = 0; v < cm->n_vc; v++)
    val_c += cm->wvc[v] * array[cm->v_ids[v]];

  return val_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct a scalar value at the cell center.
 *        array is scanned thanks to the c2v connectivity. Pointer is already
 *        located at the beginning of the cell sequence, i.e. a shift equal to
 *        c2v->idx[cm->c_id] has been done.
 *
 * \param[in] cm        pointer to a cs_cell_mesh_t structure
 * \param[in] array     array of values for each couple (v,c)
 *
 * \return the reconstructed value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_vbyc2c(const cs_cell_mesh_t      *cm,
                         const cs_real_t           *array)
{
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  cs_real_t  val_c = 0;

  if (array == NULL)
    return val_c;

  for (short int v = 0; v < cm->n_vc; v++)
    val_c += cm->wvc[v] * array[v];

  return val_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct a scalar value at the cell center.
 *        array is scanned thanks to the c2e connectivity. Pointer is already
 *        located at the beginning of the cell sequence, i.e. a shift equal to
 *        c2e->idx[cm->c_id] has been done.
 *
 * \param[in] cm        pointer to a cs_cell_mesh_t structure
 * \param[in] array     array of values for each couple (v,c)
 *
 * \return the reconstructed value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_ebyc2c(const cs_cell_mesh_t      *cm,
                         const cs_real_t           *array)
{
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PEC));

  cs_real_t  val_c = 0;

  if (array == NULL)
    return val_c;

  for (short int e = 0; e < cm->n_ec; e++)
    val_c += cm->pvol_e[e] * array[e];
  val_c /= cm->vol_c;

  return val_c;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RECO_CW_H__ */
