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

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(__cplusplus)

/*----------------------------------------------------------------------------*/
/*
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \param[in, out]  a             pointer to matrix structure
 * \param[in]       f             pointer to field, or null
 * \param[in]       iconvp        indicator
 *                                 - 1 advection
 *                                 - 0 otherwise
 * \param[in]       idiffp        indicator
 *                                 - 1 diffusion
 *                                 - 0 otherwise
 * \param[in]       ndircp        number of Dirichlet BCs
 * \param[in]       thetap        time scheme parameter
 * \param[in]       relaxp        relaxation coefficient (if < 1)
 * \param[in]       imucp         1 for temperature (with Cp), 0 otherwise
 * \param[in]       bc_coeffs     boundary condition structure
 * \param[in]       rovsdt        implicit terms (rho / dt)
 * \param[in]       i_massflux    mass flux at interior faces
 * \param[in]       b_massflux    mass flux at border faces
 * \param[in]       i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                 at interior faces for the matrix
 * \param[in]       b_visc        \f$ S_\fib \f$
 *                                 at border faces for the matrix
 * \param[in]       xcpp          Cp per cell, or null
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_compute_coeffs(cs_matrix_t                 *a,
                         const cs_field_t            *f,
                         int                          iconvp,
                         int                          idiffp,
                         int                          ndircp,
                         double                       thetap,
                         double                       relaxp,
                         int                          imucpp,
                         const cs_field_bc_coeffs_t  *bc_coeffs,
                         const cs_real_t              rovsdt[],
                         const cs_real_t              i_massflux[],
                         const cs_real_t              b_massflux[],
                         const cs_real_t              i_visc[],
                         const cs_real_t              b_visc[],
                         const cs_real_t              xcpp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the diagonal of the advection/diffusion matrix
 * for determining the variable time step, flow, Fourier.
 *
 * \tparam  stride 3 for vectors, 6 for tensors
 *
 * \param[in, out]  a                   pointer to matrix structure
 * \param[in]       f                    pointer to field, or null
 * \param[in]       iconvp               indicator
 *                                         - 1 advection
 *                                         - 0 otherwise
 * \param[in]       idiffp               indicator
 *                                         - 1 diffusion
 *                                         - 0 otherwise
 * \param[in]       tensorial_diffusion  indicator
 * \param[in]       ndircp               number of Dirichlet BCs
 * \param[in]       thetap               time scheme parameter
 * \param[in]       relaxp               relaxation coefficient (if < 1)
 * \param[in]       eb_size              extra-diagonal block size
 *                                       (1 or 3 for stride 3, 1 for stride 6)
 * \param[in]       bc_coeffs            boundary conditions structure
 * \param[in]       fimp                 implicit terms, or null
 * \param[in]       i_massflux           mass flux at interior faces
 * \param[in]       b_massflux           mass flux at border faces
 * \param[in]       i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                 at interior faces for the matrix
 * \param[in]       b_visc        \f$ S_\fib \f$
 *                                 at boundary faces for the matrix
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_matrix_compute_coeffs
(
  cs_matrix_t                 *a,
  const cs_field_t            *f,
  int                          iconvp,
  int                          idiffp,
  int                          tensorial_diffusion,
  int                          ndircp,
  cs_lnum_t                    eb_size,
  double                       thetap,
  double                       relaxp,
  const cs_field_bc_coeffs_t  *bc_coeffs,
  const cs_real_t              fimp[][stride][stride],
  const cs_real_t              i_massflux[],
  const cs_real_t              b_massflux[],
  const cs_real_t              i_visc[],
  const cs_real_t              b_visc[]
);

/*----------------------------------------------------------------------------
 * Compute legacy matrix coefficients
 *----------------------------------------------------------------------------*/

void
cs_matrix_wrapper(int                         iconvp,
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
 * Compute legacy matrix coefficients
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_matrix_wrapper(int                         iconvp,
                  int                         idiffp,
                  int                         tensorial_diffusion,
                  int                         ndircp,
                  int                         isym,
                  cs_lnum_t                   eb_size,
                  double                      thetap,
                  const cs_field_bc_coeffs_t *bc_coeffs_v,
                  const cs_real_t             fimp[][stride][stride],
                  const cs_real_t             i_massflux[],
                  const cs_real_t             b_massflux[],
                  const cs_real_t             i_visc[],
                  const cs_real_t             b_visc[],
                  cs_real_t                   da[][stride][stride],
                  cs_real_t                   xa[]);

/*----------------------------------------------------------------------------*/

#endif //defined(__cplusplus)

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
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
                    cs_real_t                  *da);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_BUILDING_H__ */
