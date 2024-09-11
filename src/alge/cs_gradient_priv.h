#pragma once

/*============================================================================
 * Private functions for gradient reconstruction.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_halo_perio.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Type for symmetric least-squares covariance matrices
   as they are adimensional, single-precision should be usable here */

typedef cs_real_t  cs_cocg_t;
typedef cs_real_t  cs_cocg_6_t[6];
typedef cs_real_t  cs_cocg_33_t[3][3];

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Semi-private inline functions
 *============================================================================*/

#if defined(__cplusplus)

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------
 * Synchronize strided gradient ghost cell values on accelerator device.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   halo_type      <-- halo type (extended or not)
 *   grad           --> gradient of a variable
 *----------------------------------------------------------------------------*/

inline static void
cs_sync_scalar_gradient_halo_d(const cs_mesh_t         *m,
                               cs_halo_type_t           halo_type,
                               cs_real_t (*restrict grad)[3])
{
  if (m->halo != NULL) {
    cs_halo_sync_d(m->halo, halo_type, CS_REAL_TYPE, 3, (cs_real_t *)grad);

    if (m->have_rotation_perio) {
      cs_sync_d2h((void  *)grad);
      cs_halo_perio_sync_var_vect(m->halo, halo_type, (cs_real_t *)grad, 3);
      cs_sync_h2d((void  *)grad);
    }
  }
}

/*----------------------------------------------------------------------------
 * Synchronize strided gradient ghost cell values on accelerator device.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   halo_type      <-- halo type (extended or not)
 *   grad           --> gradient of a variable
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
cs_sync_strided_gradient_halo_d(const cs_mesh_t         *m,
                                cs_halo_type_t           halo_type,
                                cs_real_t (*restrict grad)[stride][3])
{
  if (m->halo != NULL) {
    cs_halo_sync_d(m->halo, halo_type, CS_REAL_TYPE, stride*3,
                   (cs_real_t *)grad);

    if (m->have_rotation_perio) {
      cs_sync_d2h((void  *)grad);
      if (stride == 1)
        cs_halo_perio_sync_var_vect(m->halo, halo_type, (cs_real_t *)grad, 3);
      else if (stride == 3)
        cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
      else if (stride == 6)
        cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                             halo_type,
                                             (cs_real_t *)grad);
      cs_sync_h2d((void  *)grad);
    }
  }
}

#endif /* defined(HAVE_ACCEL) */

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

#if defined(HAVE_CUDA)

/*----------------------------------------------------------------------------
 * Compute cell gradient using least-squares reconstruction for non-orthogonal
 * meshes (nswrgp > 1).
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * cocg is computed to account for variable B.C.'s (flux).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   recompute_cocg <-- flag to recompute cocg
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   bc_coeffs      <-- B.C. structure for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or NULL
 *   cocgb          <-- saved boundary cell covariance array (on device)
 *   cocg           <-> associated cell covariance array (on device)
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

void
cs_gradient_scalar_lsq_cuda(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            cs_halo_type_t                halo_type,
                            bool                          recompute_cocg,
                            cs_real_t                     inc,
                            const cs_field_bc_coeffs_t   *bc_coeffs,
                            const cs_real_t               pvar[],
                            const cs_real_t     *restrict c_weight,
                            cs_cocg_6_t         *restrict cocgb,
                            cs_cocg_6_t         *restrict cocg,
                            cs_real_3_t         *restrict grad);

/*----------------------------------------------------------------------------
 * Compute cell gradient of a vector or tensor using least-squares
 * reconstruction for non-orthogonal meshes.
 *
 * template parameters:
 *   e2n           type of assembly algorithm used
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   n_c_iter_max   <-- maximum number of iterations for boundary correction
 *   c_eps          <-- relative tolerance for boundary correction
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable, or NULL
 *   cocgb          <-- saved boundary cell covariance array (on device)
 *   cocg           <-> cocg covariance matrix for given cell
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_gradient_strided_lsq_cuda
(
 const cs_mesh_t               *m,
 const cs_mesh_adjacencies_t   *madj,
 const cs_mesh_quantities_t    *fvq,
 const cs_halo_type_t           halo_type,
 int                            inc,
 int                            n_c_iter_max,
 cs_real_t                      c_eps,
 const cs_real_t                coefav[][stride],
 const cs_real_t                coefbv[][stride][stride],
 const cs_real_t                pvar[][stride],
 const cs_real_t               *c_weight,
 const cs_cocg_6_t             *cocgb,
 cs_cocg_6_t                   *cocg,
 cs_real_t                      grad[][stride][3]
);

/*----------------------------------------------------------------------------
 * Green-Gauss reconstruction of the gradient of a vector or tensor using
 * an initial gradient of this quantity (typically lsq).
 *
 * parameters:
 *   m                 <-- pointer to associated mesh structure
 *   fvq               <-- pointer to associated finite volume quantities
 *   cpl               <-- structure associated with internal coupling, or NULL
 *   inc               <-- if 0, solve on increment; 1 otherwise
 *   porous_model      <-- type of porous model used
 *   warped_correction <-- apply warped faces correction ?
 *   coefav            <-- B.C. coefficients for boundary face normals
 *   coefbv            <-- B.C. coefficients for boundary face normals
 *   pvar              <-- variable
 *   c_weight          <-- weighted gradient coefficient variable
 *   r_grad            <-- gradient used for reconstruction
 *   grad              --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_gradient_strided_gg_r_cuda
(
 const cs_mesh_t              *m,
 const cs_mesh_adjacencies_t  *madj,
 const cs_mesh_quantities_t   *fvq,
 cs_halo_type_t                halo_type,
 int                           inc,
 int                           porous_model,
 bool                          warped_correction,
 const cs_real_t               coefav[][stride],
 const cs_real_t               coefbv[][stride][stride],
 const cs_real_t               pvar[][stride],
 const cs_real_t              *c_weight,
 const cs_real_t               r_grad[][stride][3],
 cs_real_t                     grad[][stride][3]
);

#endif /* defined(HAVE_CUDA) */

#endif /* defined(__cplusplus) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
