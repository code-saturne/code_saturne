#ifndef __CS_FINITE_VOLUME_OPERATORS_H__
#define __CS_FINITE_VOLUME_OPERATORS_H__

/*============================================================================
 * Finite volume operators
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_divmas
 *----------------------------------------------------------------------------*/

void CS_PROCF (divmas, DIVMAS)
(
 const cs_int_t  *const   n_cells_ext,
 const cs_int_t  *const   n_cells,
 const cs_int_t  *const   n_i_faces,
 const cs_int_t  *const   n_b_faces,
 const cs_int_t  *const   init,
 const cs_int_t  *const   nfecra,
 const cs_lnum_2_t        i_face_cells[],
 const cs_int_t           b_face_cells[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 cs_real_t                diverg[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_itrmas
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmas, ITRMAS)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_itrgrp
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrp, ITRGRP)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
 const cs_int_t  *const   nfecra,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 cs_real_3_t              frcxt[],
 cs_real_t                pvar[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[],
 cs_real_t                diverg[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_projts
 *----------------------------------------------------------------------------*/

void CS_PROCF (projts, PROJTS)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_int_t  *const   nfecra,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[]);

/*----------------------------------------------------------------------------
 * Wrapper to cs_finite_volume_projtv
 *----------------------------------------------------------------------------*/

void CS_PROCF (projtv, PROJTV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   nfecra,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_6_t        viscel[],
 const cs_real_2_t        weighf[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*! \brief This function adds the integrated mass flux on the cells.

  \f[
  \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
  \f]

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     n_cells_ext   number of extended (real + ghost) cells
 * \param[in]     n_cells       number of cells
 * \param[in]     n_i_faces     number of interior faces
 * \param[in]     n_b_faces     number of boundary faces
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_face_cells  cell indexes of interior faces
 * \param[in]     b_face_cells  cell indexes of boundary faces
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_divergence(
                            int                       n_cells_ext,
                            int                       n_cells,
                            int                       n_i_faces,
                            int                       n_b_faces,
                            int                       init,
                            const cs_lnum_2_t         i_face_cells[],
                            const cs_int_t            b_face_cells[],
                            const cs_real_t           i_massflux[],
                            const cs_real_t           b_massflux[],
                            cs_real_t       *restrict diverg);

/*----------------------------------------------------------------------------*/

/*! \brief This function updates the face mass flux with the face pressure
  (or pressure increment, or pressure double increment) gradient.

   \f[
   \dot{m}_\ij = \dot{m}_\ij
               - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
   \f]

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_face_gradient_scalar(
                                      int                       init,
                                      int                       inc,
                                      int                       imrgra,
                                      int                       iccocg,
                                      int                       nswrgp,
                                      int                       imligp,
                                      int                       iphydp,
                                      int                       iwarnp,
                                      double                    epsrgp,
                                      double                    climgp,
                                      double                    extrap,
                                      cs_real_3_t     *restrict frcxt,
                                      cs_real_t       *restrict pvar,
                                      const cs_real_t           coefap[],
                                      const cs_real_t           coefbp[],
                                      const cs_real_t           cofafp[],
                                      const cs_real_t           cofbfp[],
                                      const cs_real_t           i_visc[],
                                      const cs_real_t           b_visc[],
                                      const cs_real_t           viselx[],
                                      const cs_real_t           visely[],
                                      const cs_real_t           viselz[],
                                      cs_real_t       *restrict i_massflux,
                                      cs_real_t       *restrict b_massflux);

/*----------------------------------------------------------------------------*/

/*! \brief This function updates the cell mass flux divergence with the face
  pressure (or pressure increment, or pressure double increment) gradient.

  \f[
  \dot{m}_\ij = \dot{m}_\ij
              - \sum_j \Delta t \grad_\fij p \cdot \vect{S}_\ij
  \f]

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                 (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_div_face_gradient_scalar(
                                          int                       init,
                                          int                       inc,
                                          int                       imrgra,
                                          int                       iccocg,
                                          int                       nswrgp,
                                          int                       imligp,
                                          int                       iphydp,
                                          int                       iwarnp,
                                          double                    epsrgp,
                                          double                    climgp,
                                          double                    extrap,
                                          cs_real_3_t     *restrict frcxt,
                                          cs_real_t       *restrict pvar,
                                          const cs_real_t           coefap[],
                                          const cs_real_t           coefbp[],
                                          const cs_real_t           cofafp[],
                                          const cs_real_t           cofbfp[],
                                          const cs_real_t           i_visc[],
                                          const cs_real_t           b_visc[],
                                          const cs_real_t           viselx[],
                                          const cs_real_t           visely[],
                                          const cs_real_t           viselz[],
                                          cs_real_t       *restrict diverg);

/*----------------------------------------------------------------------------*/

/*! \brief This function projects the external source terms to the faces
  in coherence with itrmas.f90 for the improved hydrostatic pressure
  algorithm (iphydr=1).

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 _____________________________________________________________________________*/
/*!
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgu        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_face_source_terms_scalar(
                                          int                       init,
                                          int                       nswrgu,
                                          const cs_real_3_t         frcxt[],
                                          const cs_real_t           cofbfp[],
                                          cs_real_t       *restrict i_massflux,
                                          cs_real_t       *restrict b_massflux,
                                          const cs_real_t           i_visc[],
                                          const cs_real_t           b_visc[],
                                          const cs_real_t           viselx[],
                                          const cs_real_t           visely[],
                                          const cs_real_t           viselz[]);

/*----------------------------------------------------------------------------*/

/*! \brief This function projects the external source terms to the faces
  in coherence with itrmav.f90 for the improved hydrostatic pressure
  algorithm (iphydr=1).

*/
/*------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 ______________________________________________________________________________!
 * \param[in]     init           indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_finite_volume_face_source_terms_vector(
                                          int                       init,
                                          int                       nswrgp,
                                          int                       ircflp,
                                          const cs_real_3_t         frcxt[],
                                          const cs_real_t           cofbfp[],
                                          const cs_real_t           i_visc[],
                                          const cs_real_t           b_visc[],
                                          const cs_real_6_t         viscel[],
                                          const cs_real_2_t         weighf[],
                                          cs_real_t       *restrict i_massflux,
                                          cs_real_t       *restrict b_massflux);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FINITE_VOLUME_OPERATORS_H__ */
