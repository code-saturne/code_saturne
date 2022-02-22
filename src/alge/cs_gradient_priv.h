#ifndef __CS_GRADIENT_CUDA_H__
#define __CS_GRADIENT_CUDA_H__

/*============================================================================
 * Private functions for gradient reconstruction.
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
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   fext           <-- exterior force generating pressure
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
                            int                           hyd_p_flag,
                            cs_real_t                     inc,
                            const cs_real_3_t             f_ext[],
                            const cs_real_t               coefap[],
                            const cs_real_t               coefbp[],
                            const cs_real_t               pvar[],
                            const cs_real_t     *restrict c_weight,
                            cs_cocg_6_t         *restrict cocg,
                            cs_cocg_6_t         *restrict cocgb,
                            cs_real_3_t         *restrict grad);

#endif /* defined(HAVE_CUDA) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GRADIENT_CUDA_H__ */
