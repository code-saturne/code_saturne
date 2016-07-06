#ifndef __CS_CDO_SCHEME_GEOMETRY_H__
#define __CS_CDO_SCHEME_GEOMETRY_H__

/*============================================================================
 * Geometric computations for building discretization operators which is
 * shared by several files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"

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
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_get_face_wbs0(short int                   f,
                     const cs_cdo_locmesh_t     *lm,
                     cs_real_t                  *wvf,
                     cs_real_t                  *pefc_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 *
 * \return the volume of p_{f,c}
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_get_face_wbs1(short int                   f,
                     const cs_cdo_locmesh_t     *lm,
                     cs_real_t                  *wvf,
                     cs_real_t                  *pefc_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] grd_c      gradient of the Lagrange function related to xc
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_get_face_wbs2(short int                   f,
                     const cs_cdo_locmesh_t     *lm,
                     cs_real_3_t                 grd_c,
                     cs_real_t                  *wvf,
                     cs_real_t                  *pefc_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] grd_c      gradient of the Lagrange function related to xc
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 *
 * \return the volume of p_{f,c}
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_get_face_wbs3(short int                   f,
                     const cs_cdo_locmesh_t     *lm,
                     cs_real_3_t                 grd_c,
                     cs_real_t                  *wvf,
                     cs_real_t                  *pefc_vol);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_SCHEME_GEOMETRY_H__ */
