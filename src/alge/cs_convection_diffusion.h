#ifndef CS_CONVECTION_DIFFUSION_H
#define CS_CONVECTION_DIFFUSION_H

/*============================================================================
 * Convection-diffusion operators.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "base/cs_halo.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh_quantities.h"
#include "cdo/cs_equation_param.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * NVD/TVD Advection Scheme
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_NVD_GAMMA      = 0,      /* GAMMA            */
  CS_NVD_SMART      = 1,      /* SMART            */
  CS_NVD_CUBISTA    = 2,      /* CUBISTA          */
  CS_NVD_SUPERBEE   = 3,      /* SUPERBEE         */
  CS_NVD_MUSCL      = 4,      /* MUSCL            */
  CS_NVD_MINMOD     = 5,      /* MINMOD           */
  CS_NVD_CLAM       = 6,      /* CLAM             */
  CS_NVD_STOIC      = 7,      /* STOIC            */
  CS_NVD_OSHER      = 8,      /* OSHER            */
  CS_NVD_WASEB      = 9,      /* WASEB            */
  CS_NVD_VOF_HRIC   = 10,     /* M-HRIC for VOF   */
  CS_NVD_VOF_CICSAM = 11,     /* M-CICSAM for VOF */
  CS_NVD_VOF_STACS  = 12,     /* STACS for VOF    */
  CS_NVD_N_TYPES    = 13      /* number of NVD schemes */

} cs_nvd_type_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointer to slope test indicator field values if active.
 *----------------------------------------------------------------------------*/

cs_real_t *
cs_get_v_slope_test(int                        f_id,
                    const cs_equation_param_t  eqp);

/*----------------------------------------------------------------------------
 * Compute the beta blending coefficient of the
 * beta limiter (ensuring preservation of a given min/max pair of values).
 *----------------------------------------------------------------------------*/

void
cs_beta_limiter_building(int              f_id,
                         int              inc,
                         const cs_real_t  rovsdt[]);

/*----------------------------------------------------------------------------
 * Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field.
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_scalar(const cs_field_t           *f,
                               const cs_equation_param_t   eqp,
                               int                         icvflb,
                               int                         inc,
                               int                         imasac,
                               const cs_real_t            *pvar,
                               const int                   icvfli[],
                               cs_field_bc_coeffs_t       *bc_coeffs,
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               const cs_real_t            *c_weight,
                               cs_real_t                  *rhs,
                               cs_real_2_t                 i_flux[],
                               cs_real_t                   b_flux[]);

/*----------------------------------------------------------------------------*
 * Add the explicit part of the convection/diffusion terms of a transport
 * equation of a scalar field such as the temperature.
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_thermal(const cs_field_t           *f,
                                const cs_equation_param_t   eqp,
                                int                         inc,
                                int                         imasac,
                                const cs_real_t            *pvar,
                                cs_field_bc_coeffs_t       *bc_coeffs,
                                const cs_real_t             i_massflux[],
                                const cs_real_t             b_massflux[],
                                const cs_real_t             i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t            *c_weight,
                                const cs_real_t             xcpp[],
                                cs_real_t                  *rhs);

/*----------------------------------------------------------------------------
 * Update face flux with convection contribution of a standard transport
 *----------------------------------------------------------------------------*/

void
cs_face_convection_scalar(int                         idtvar,
                          int                         f_id,
                          const cs_equation_param_t   eqp,
                          int                         icvflb,
                          int                         inc,
                          cs_real_t                  *pvar,
                          const cs_real_t            *pvara,
                          const int                   icvfli[],
                          cs_field_bc_coeffs_t       *bc_coeffs,
                          const cs_real_t             i_massflux[],
                          const cs_real_t             b_massflux[],
                          cs_real_2_t                 i_conv_flux[],
                          cs_real_t                   b_conv_flux[]);

/*----------------------------------------------------------------------------
 * Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field.
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_scalar_at_faces
(
  const cs_field_t           *f,
  const cs_equation_param_t   eqp,
  int                         icvflb,
  cs_lnum_t                   n_i_faces,
  cs_lnum_t                   n_b_faces,
  const cs_lnum_t            *i_face_ids,
  const cs_lnum_t            *b_face_ids,
  const cs_real_t            *restrict pvar,
  const int                   icvfli[],
  cs_field_bc_coeffs_t       *bc_coeffs,
  const cs_real_t             i_massflux[],
  const cs_real_t             b_massflux[],
  const cs_real_t             i_visc[],
  const cs_real_t             b_visc[],
  const cs_real_t            *c_weight,
  const cs_real_t            *xcpp,
  cs_real_2_t                 i_flux[],
  cs_real_t                   b_flux[]
);

/*----------------------------------------------------------------------------
 * Add the explicit part of the convection/diffusion terms of a transport
 * equation of a vector field \f$ \vect{\varia} \f$.
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_vector(int                         idtvar,
                               int                         f_id,
                               const cs_equation_param_t   eqp,
                               int                         icvflb,
                               int                         inc,
                               int                         ivisep,
                               int                         imasac,
                               cs_real_3_t                *pvar,
                               const cs_real_3_t          *pvara,
                               const int                   icvfli[],
                               cs_field_bc_coeffs_t       *bc_coeffs,
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               const cs_real_t             i_secvis[],
                               const cs_real_t             b_secvis[],
                               cs_real_3_t                *i_pvar,
                               cs_real_3_t                *b_pvar,
                               cs_real_3_t                *rhs);

/*----------------------------------------------------------------------------
 * Add the explicit part of the convection/diffusion terms of a transport
 * equation of a tensor field \f$ \tens{\varia} \f$.
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_tensor(int                          idtvar,
                               int                          f_id,
                               const cs_equation_param_t    eqp,
                               int                          icvflb,
                               int                          inc,
                               int                          imasac,
                               cs_real_6_t                 *pvar,
                               const cs_real_6_t           *pvara,
                               cs_field_bc_coeffs_t        *bc_coeffs,
                               const cs_real_t              i_massflux[],
                               const cs_real_t              b_massflux[],
                               const cs_real_t              i_visc[],
                               const cs_real_t              b_visc[],
                               cs_real_6_t                 *rhs);

/*----------------------------------------------------------------------------
 * Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_scalar(int                         idtvar,
                                int                         f_id,
                                const cs_equation_param_t   eqp,
                                int                         inc,
                                const cs_real_t            *pvar,
                                const cs_real_t            *pvara,
                                cs_field_bc_coeffs_t       *bc_coeffs,
                                const cs_real_t             i_visc[],
                                const cs_real_t             b_visc[],
                                cs_real_6_t                *viscel,
                                const cs_real_2_t           weighf[],
                                const cs_real_t             weighb[],
                                cs_real_t                  *rhs);

/*-----------------------------------------------------------------------------
 * Add explicit part of the terms of diffusion by a left-multiplying
 * symmetric tensorial diffusivity for a transport equation of a vector field
 * \f$ \vect{\varia} \f$.
 *----------------------------------------------------------------------------*/

void
cs_anisotropic_left_diffusion_vector(int                         idtvar,
                                     int                         f_id,
                                     const cs_equation_param_t   eqp,
                                     int                         inc,
                                     int                         ivisep,
                                     cs_real_3_t                *pvar,
                                     const cs_real_3_t          *pvara,
                                     cs_field_bc_coeffs_t       *bc_coeffs,
                                     const cs_real_33_t          i_visc[],
                                     const cs_real_t             b_visc[],
                                     const cs_real_t             i_secvis[],
                                     cs_real_3_t                *rhs);

/*-----------------------------------------------------------------------------
 * Add explicit part of the terms of diffusion by a right-multiplying
 * symmetric tensorial diffusivity for a transport equation of a vector field
 * \f$ \vect{\varia} \f$.
 *----------------------------------------------------------------------------*/

void
cs_anisotropic_right_diffusion_vector(int                          idtvar,
                                      int                          f_id,
                                      const cs_equation_param_t    eqp,
                                      int                          inc,
                                      cs_real_3_t                 *pvar,
                                      const cs_real_3_t           *pvara,
                                      cs_field_bc_coeffs_t        *bc_coeffs,
                                      const cs_real_t              i_visc[],
                                      const cs_real_t              b_visc[],
                                      cs_real_6_t                 *viscel,
                                      const cs_real_2_t            weighf[],
                                      const cs_real_t              weighb[],
                                      cs_real_3_t                 *rhs);

/*----------------------------------------------------------------------------
 * Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_tensor(int                          idtvar,
                                int                          f_id,
                                const cs_equation_param_t    eqp,
                                int                          inc,
                                cs_real_6_t                 *pvar,
                                const cs_real_6_t           *pvara,
                                cs_field_bc_coeffs_t        *bc_coeffs,
                                const cs_real_t              i_visc[],
                                const cs_real_t              b_visc[],
                                cs_real_6_t                 *viscel,
                                const cs_real_2_t            weighf[],
                                const cs_real_t              weighb[],
                                cs_real_6_t                 *rhs);

/*----------------------------------------------------------------------------
 * Update the face mass flux with the face pressure (or pressure
 * increment, or pressure double increment) gradient.
 *----------------------------------------------------------------------------*/

void
cs_face_diffusion_potential(const cs_field_t           *f,
                            const cs_equation_param_t  *eqp,
                            const cs_mesh_t            *m,
                            cs_mesh_quantities_t       *fvq,
                            int                         init,
                            int                         inc,
                            int                         iphydp,
                            cs_real_3_t                *frcxt,
                            cs_real_t                  *pvar,
                            cs_field_bc_coeffs_t       *bc_coeffs,
                            const cs_real_t             i_visc[],
                            const cs_real_t             b_visc[],
                            cs_real_t                  *visel,
                            cs_real_t                  *i_massflux,
                            cs_real_t                  *b_massflux);

/*----------------------------------------------------------------------------
 * Add the explicit part of the pressure gradient term to the mass flux
 * in case of anisotropic diffusion of the pressure field \f$ P \f$.
 *----------------------------------------------------------------------------*/

void
cs_face_anisotropic_diffusion_potential
(
  const cs_field_t           *f,
  const cs_equation_param_t  *eqp,
  const cs_mesh_t            *m,
  cs_mesh_quantities_t       *fvq,
  int                         init,
  int                         inc,
  int                         iphydp,
  cs_real_3_t                *frcxt,
  cs_real_t                  *pvar,
  cs_field_bc_coeffs_t       *bc_coeffs,
  const cs_real_t             i_visc[],
  const cs_real_t             b_visc[],
  cs_real_6_t                *viscel,
  const cs_real_2_t           weighf[],
  const cs_real_t             weighb[],
  cs_real_t                  *i_massflux,
  cs_real_t                  *b_massflux
);

/*----------------------------------------------------------------------------
 * Update the cell mass flux divergence with the face pressure (or
 * pressure increment, or pressure double increment) gradient.
 *----------------------------------------------------------------------------*/

void
cs_diffusion_potential(const cs_field_t           *f,
                       const cs_equation_param_t  *eqp,
                       const cs_mesh_t            *m,
                       cs_mesh_quantities_t       *fvq,
                       int                         init,
                       int                         inc,
                       int                         iphydp,
                       cs_real_3_t                *frcxt,
                       cs_real_t                  *pvar,
                       cs_field_bc_coeffs_t       *bc_coeffs,
                       const cs_real_t             i_visc[],
                       const cs_real_t             b_visc[],
                       cs_real_t                   visel[],
                       cs_real_t                  *diverg);

/*----------------------------------------------------------------------------
 * Add the explicit part of the divergence of the mass flux due to the
 * pressure gradient (analog to cs_anisotropic_diffusion_scalar).
 *----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_potential(const cs_field_t           *f,
                                   const cs_equation_param_t  *eqp,
                                   const cs_mesh_t            *m,
                                   cs_mesh_quantities_t       *fvq,
                                   int                         init,
                                   int                         inc,
                                   int                         iphydp,
                                   cs_real_3_t       *restrict frcxt,
                                   cs_real_t         *restrict pvar,
                                   cs_field_bc_coeffs_t       *bc_coeffs,
                                   const cs_real_t             i_visc[],
                                   const cs_real_t             b_visc[],
                                   cs_real_6_t       *restrict viscel,
                                   const cs_real_2_t           weighf[],
                                   const cs_real_t             weighb[],
                                   cs_real_t         *restrict diverg);

/*----------------------------------------------------------------------------
 * Compute the upwind gradient used in the slope tests.
 *----------------------------------------------------------------------------*/

template <typename T>
void
cs_slope_test_gradient(int                         f_id,
                       cs_dispatch_context        &ctx,
                       const T                   (*grad)[3],
                       T                         (*grdpa)[3],
                       const cs_real_t            *pvar,
                       const cs_real_t             val_f[],
                       const cs_real_t            *i_massflux);

/*----------------------------------------------------------------------------
 * Compute the upwind gradient used in the slope tests.
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride, typename T>
void
cs_slope_test_gradient_strided
  (cs_dispatch_context         &ctx,
   const T                      grad[][stride][3],
   T                 (*restrict grdpa)[stride][3],
   const cs_real_t              pvar[][stride],
   const cs_real_t              val_f[][stride],
   const cs_real_t             *i_massflux);

/*----------------------------------------------------------------------------
 * Compute the upwind gradient used in the pure SOLU schemes
 * (observed in the litterature).
 *----------------------------------------------------------------------------*/

template <typename T>
void
cs_upwind_gradient(cs_dispatch_context          &ctx,
                   const int                     inc,
                   const cs_halo_type_t          halo_type,
                   const cs_field_bc_coeffs_t   *bc_coeffs,
                   const cs_real_t               i_massflux[],
                   const cs_real_t               b_massflux[],
                   const cs_real_t              *pvar,
                   T                           (*grdpa)[3]);

/*----------------------------------------------------------------------------
 * Compute the upwind gradient used in the pure SOLU schemes
 * (observed in the litterature) for a vector or tensor
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride, typename T>
void
cs_upwind_gradient_strided(cs_dispatch_context          &ctx,
                           const int                     inc,
                           const cs_halo_type_t          halo_type,
                           const cs_field_bc_coeffs_t   *bc_coeffs,
                           const cs_real_t               i_massflux[],
                           const cs_real_t               b_massflux[],
                           const cs_real_t     *restrict pvar[stride],
                           T                  (*restrict grdpa)[stride][3]);

/*----------------------------------------------------------------------------
 * Compute the local cell Courant number as the maximum of all cell face based
 * Courant number at each cell.
 *----------------------------------------------------------------------------*/

void
cs_cell_courant_number(const cs_field_t     *f,
                       cs_dispatch_context  &ctx,
                       cs_real_t            *courant);

/*----------------------------------------------------------------------------
 * Compute balance contribution of the transpose grad(vel) term
 * and grad(-2/3 div(vel))
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_secvis
(
  cs_dispatch_context         &ctx,
  const cs_mesh_t             *m,
  const cs_mesh_quantities_t  *mq,
  cs_real_t                    thetap,
  const cs_real_t              i_visc[],
  const cs_real_t              i_secvis[],
  const cs_real_t              b_secvis[],
  const cs_rreal_t             gradv[][3][3],
  cs_real_3_t        *restrict rhs
);

/*----------------------------------------------------------------------------
 * Compute balance contribution of the transpose grad(vel) term
 * and grad(-2/3 div(vel)) with anisotropic.
 *----------------------------------------------------------------------------*/

void
cs_convection_anisotropic_leff_diffusion_secvis
(
  cs_dispatch_context         &ctx,
  const cs_mesh_t             *m,
  const cs_mesh_quantities_t  *mq,
  const cs_real_33_t           i_visc[],
  const cs_real_t              i_secvis[],
  const cs_real_t              gradv[][3][3],
  cs_real_3_t        *restrict rhs
);

/*----------------------------------------------------------------------------
 * Query convection-diffusion scheme variants.
 *----------------------------------------------------------------------------*/

int
cs_convection_diffusion_get_scheme_version(void);

/*----------------------------------------------------------------------------
 * Allow reverting to older convection-diffusion scheme variants.
 *----------------------------------------------------------------------------*/

void
cs_convection_diffusion_set_scheme_version(int  version);

/*----------------------------------------------------------------------------*/

#endif /* CS_CONVECTION_DIFFUSION_H */
