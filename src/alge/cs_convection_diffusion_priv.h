#ifndef __CS_CONVECTION_DIFFUSION_CUDA_H__
#define __CS_CONVECTION_DIFFUSION_CUDA_H__

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


/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

#if defined(HAVE_CUDA)

void
cs_convection_diffusion_vector_cuda(const cs_mesh_t             *mesh,
                                    const cs_mesh_quantities_t  *fvq,
                                    const cs_real_3_t  *restrict pvar,
                                    const cs_real_t              i_massflux[],
                                    const cs_real_33_t          *grad,
                                    cs_real_33_t                *grdpa,
                                    cs_real_3_t       *restrict  rhs,
                                    const cs_real_3_t  *restrict coefav,
                                    const cs_real_33_t *restrict coefbv,
                                    const int                    inc,
                                    const bool                   flag1,
                                    const bool                   flag2,
                                    const bool                   perf);

#endif

/* defined(HAVE_CUDA) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

END_C_DECLS
#endif /* __CS_GRADIENT_CUDA_H__ */
