#ifndef __CS_GRADIENT_BOUNDARY_H__
#define __CS_GRADIENT_BOUNDARY_H__

/*============================================================================
 * Gradient reconstruction at boundaries.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_internal_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a scalar at boundary face I' positions
 *         using least-squares interpolation.
 *
 * This assumes ghost cell values for the variable (var) are up-to-date.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neumann BC for a given face.
 *
 * We still use the provided BC's when possible, for the following cases:
 * - Given a non-uniform Dirichlet condition and a non-orthogonal mesh,
 *   the Dirichlet values at face centers (shifted by II' relative to I')
 *   can convey a portion of the information of the gradient along II'.
 * - For cells with multiple boundary faces, information from faces whose
 *   normals are not orthogonal to II' can also provide a significant
 *   contribution to the normal.
 *
 * For coupled faces (with internal coupling), we assume a Neumann BC,
 * as the values that could be computed at the intersections of the
 * segments joining coupled cell centers and the matching boundary face
 * do not currently account for the non-linearity associated with wall laws.
 *
 * A more precise computation for cells with multiple coupled boundary faces
 * would be possible by using multiple passes, first estimating boundary face
 * values without recontruction, then using that information for faces which
 * are not tangential to II' in a second pass. Note though that as the mesh
 * is refined, the portion of the boundary adjacent to cells with
 * multiple coupled faces (i.e. at corners and ridges) diminishes.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_boundary_iprime_lsq_s(const cs_mesh_t               *m,
                                  const cs_mesh_quantities_t    *fvq,
                                  const cs_internal_coupling_t  *cpl,
                                  cs_lnum_t                      n_faces,
                                  const cs_lnum_t               *face_ids,
                                  cs_halo_type_t                 halo_type,
                                  double                         clip_coeff,
                                  const cs_real_t               *bc_coeff_a,
                                  const cs_real_t               *bc_coeff_b,
                                  const cs_real_t                c_weight[],
                                  const cs_real_t                var[],
                                  cs_real_t           *restrict  var_iprime);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a vector at boundary face I' positions
 *         using least-squares interpolation.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * This function uses a local iterative approach to compute the cell gradient,
 * as handling of the boundary condition terms b in higher dimensions
 * would otherwise require solving higher-dimensional systems, often at
 * a higher cost.
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neuman BC.
 *
 * The same logic is applied as for \ref cs_gradient_boundary_iprime_lsq_s.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_boundary_iprime_lsq_v(const cs_mesh_t               *m,
                                  const cs_mesh_quantities_t    *fvq,
                                  const cs_internal_coupling_t  *cpl,
                                  cs_lnum_t                      n_faces,
                                  const cs_lnum_t               *face_ids,
                                  cs_halo_type_t                 halo_type,
                                  double                         clip_coeff,
                                  const cs_real_t               *bc_coeff_a[3],
                                  const cs_real_t               *bc_coeff_b[3][3],
                                  const cs_real_t                c_weight[],
                                  const cs_real_t                var[][6],
                                  cs_real_t            *restrict var_iprime[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a symmetric tensor at boundary face
 *         I' positions using least-squares interpolation.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * This function uses a local iterative approach to compute the cell gradient,
 * as handling of the boundary condition terms b in higher dimensions
 * would otherwise require solving higher-dimensional systems, often at
 * a higher cost.
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neuman BC.
 *
 * The same logic is applied as for \ref cs_gradient_boundary_iprime_lsq_s.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_boundary_iprime_lsq_t(const cs_mesh_t               *m,
                                  const cs_mesh_quantities_t    *fvq,
                                  const cs_internal_coupling_t  *cpl,
                                  cs_lnum_t                      n_faces,
                                  const cs_lnum_t               *face_ids,
                                  cs_halo_type_t                 halo_type,
                                  double                         clip_coeff,
                                  const cs_real_t               *bc_coeff_a[6],
                                  const cs_real_t               *bc_coeff_b[6][6],
                                  const cs_real_t                c_weight[],
                                  const cs_real_t                var[][6],
                                  cs_real_t            *restrict var_iprime[6]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GRADIENT_BOUNDARY__ */
