#ifndef __CS_MATRIX_BUILDING_H__
#define __CS_MATRIX_BUILDING_H__

/*============================================================================
 * Matrix building
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
#include "cs_math.h"
#include "cs_parameters.h"  // for BC types

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"

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
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_scalar (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_scalar(int                         iconvp,
                         int                         idiffp,
                         int                         ndircp,
                         int                         isym,
                         double                      thetap,
                         int                         imucpp,
                         const cs_field_bc_coeffs_t *bc_coeffs,
                         const cs_real_t             rovsdt[],
                         const cs_real_t             i_massflux[],
                         const cs_real_t             b_massflux[],
                         const cs_real_t             i_visc[],
                         const cs_real_t             b_visc[],
                         const cs_real_t             xcpp[],
                         cs_real_t                   da[],
                         cs_real_t                   xa[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_vector (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_vector(int                         iconvp,
                         int                         idiffp,
                         int                         tensorial_diffusion,
                         int                         ndircp,
                         int                         isym,
                         cs_lnum_t                   eb_size,
                         double                      thetap,
                         const cs_field_bc_coeffs_t *bc_coeffs_v,
                         const cs_real_t             fimp[][3][3],
                         const cs_real_t             i_massflux[],
                         const cs_real_t             b_massflux[],
                         const cs_real_t             i_visc[],
                         const cs_real_t             b_visc[],
                         cs_real_t                   da[][3][3],
                         cs_real_t                   xa[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_matrix_tensor (or its counterpart for
 * symmetric matrices)
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper_tensor(int                         iconvp,
                         int                         idiffp,
                         int                         tensorial_diffusion,
                         int                         ndircp,
                         int                         isym,
                         double                      thetap,
                         const cs_field_bc_coeffs_t *bc_coeffs_ts,
                         const cs_real_66_t          fimp[],
                         const cs_real_t             i_massflux[],
                         const cs_real_t             b_massflux[],
                         const cs_real_t             i_visc[],
                         const cs_real_t             b_visc[],
                         cs_real_66_t                da[],
                         cs_real_t                   xa[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \param[in]     m             pointer to mesh structure
 * \param[in]     iconvp        indicator
 *                               - 1 advection
 *                               - 0 otherwise
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion
 *                               - 0 otherwise
 * \param[in]     isym          indicator
 *                               - 1 symmetric matrix
 *                               - 2 non symmmetric matrix
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at border faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_visc        \f$ S_\fib \f$
 *                               at border faces for the matrix
 * \param[out]    da            diagonal part of the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_time_step(const cs_mesh_t            *m,
                    int                         iconvp,
                    int                         idiffp,
                    int                         isym,
                    const cs_field_bc_coeffs_t *bc_coeffs,
                    const cs_real_t             i_massflux[],
                    const cs_real_t             b_massflux[],
                    const cs_real_t             i_visc[],
                    const cs_real_t             b_visc[],
                    cs_real_t         *restrict da);

END_C_DECLS

#endif /* __CS_MATRIX_BUILDING_H__ */
