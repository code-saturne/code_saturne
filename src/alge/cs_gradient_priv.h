#ifndef __CS_GRADIENT_CUDA_H__
#define __CS_GRADIENT_CUDA_H__

/*============================================================================
 * Private functions for gradient reconstruction.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include "cs_halo.h"
#include "cs_internal_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 *  Global variables
 *============================================================================*/

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
 *   coefap         <-- B.C. coefficients for boundary face normals
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or NULL
 *   cocg           <-> associated cell covariance array (on device)
 *   cocgb          <-> saved boundary cell covariance array (on device)
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

void
cs_gradient_scalar_lsq_cuda(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            cs_halo_type_t                halo_type,
                            bool                          recompute_cocg,
                            cs_real_t                     inc,
                            const cs_real_t               coefap[],
                            const cs_real_t               coefbp[],
                            const cs_real_t               pvar[],
                            const cs_real_t     *restrict c_weight,
                            cs_cocg_6_t         *restrict cocg,
                            cs_cocg_6_t         *restrict cocgb,
                            cs_real_3_t         *restrict grad);

void
cs_lsq_vector_gradient_cuda(const cs_mesh_t        *m,
                     const cs_mesh_adjacencies_t   *madj,
                     const cs_mesh_quantities_t    *fvq,
                     const cs_halo_type_t           halo_type,
                     const int                      inc,
                     const cs_real_3_t    *restrict coefav,
                     const cs_real_33_t   *restrict coefbv,
                     const cs_real_3_t    *restrict pvar,
                     const cs_real_t      *restrict c_weight,
                     cs_cocg_6_t          *restrict cocg,
                     cs_cocg_6_t          *restrict cocgb,
                     cs_real_33_t         *restrict gradv,
                     cs_real_33_t         *restrict rhs);

void
cs_reconstruct_vector_gradient_cuda(const cs_mesh_t              *m,
                                    const cs_mesh_adjacencies_t  *madj,
                                    const cs_mesh_quantities_t   *fvq,
                                    const cs_internal_coupling_t *cpl,
                                    cs_halo_type_t                halo_type,
                                    int                           inc,
                                    const cs_real_3_t   *restrict coefav,
                                    const cs_real_33_t  *restrict coefbv,
                                    const cs_real_3_t   *restrict pvar,
                                    const cs_real_t     *restrict c_weight,
                                    const cs_real_33_t        *restrict r_grad,
                                    cs_real_33_t        *restrict grad,
                                    const bool                   *coupled_faces,
                                    cs_lnum_t                     cpl_stride,
                                    bool                          test_bool,
                                    bool                          perf);

void
_gradient_vector_cuda(const cs_mesh_t     *mesh,
                      cs_real_3_t         *_bc_coeff_a,
                      cs_real_33_t        *_bc_coeff_b,
                      bool                 a_null,
                      bool                 b_null,
                      bool                 perf);
                      
#endif

/* defined(HAVE_CUDA) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

END_C_DECLS
#ifdef __cplusplus
/**
 * This template will be instantited with stride = 1, 3, 6, 9
*/
template <cs_lnum_t stride>
void
cs_lsq_vector_gradient_strided_cuda(const cs_mesh_t               *m,
                     const cs_mesh_adjacencies_t   *madj,
                     const cs_mesh_quantities_t    *fvq,
                     const cs_halo_type_t           halo_type,
                     const int                      inc,
                     const cs_real_t (*restrict coefav)[stride],
                     const cs_real_t (*restrict coefbv)[stride][stride],
                     const cs_real_t (*restrict pvar)[stride],
                     const cs_real_t      *restrict c_weight,
                     cs_cocg_6_t          *restrict cocg,
                     cs_cocg_6_t          *restrict cocgb,
                     cs_real_t (*restrict gradv)[stride][3],
                     cs_real_t (*restrict rhs)[stride][3],
                     cs_lnum_t n_c_iter_max,
                     cs_real_t c_eps);
#endif
#endif /* __CS_GRADIENT_CUDA_H__ */
