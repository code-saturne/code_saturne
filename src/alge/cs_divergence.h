#ifndef __CS_DIVERGENCE_H__
#define __CS_DIVERGENCE_H__

/*============================================================================
 * Divergence operators.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add \f$ \rho \vect{u} \cdot \vect{s}_\ij\f$ to
 * the mass flux \f$ \dot{m}_\ij \f$.
 *
 * For the reconstruction, \f$ \gradt \left(\rho \vect{u} \right) \f$ is
 * computed with the following approximated boundary conditions:
 *  - \f$ \vect{a}_{\rho u} = \rho_\fib \vect{a}_u \f$
 *  - \f$ \tens{b}_{\rho u} = \tens{b}_u \f$
 *
 * For the mass flux at the boundary we have:
 * \f[
 * \dot{m}_\ib = \left[ \rho_\fib \vect{a}_u  + \rho_\fib \tens{b}_u \vect{u}
 * + \tens{b}_u \left(\gradt \vect{u} \cdot \vect{\centi \centip}\right)\right]
 * \cdot \vect{s}_\ij
 * \f]
 * The last equation uses some approximations detailed in the theory guide.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     f_id          field id (or -1)
 * \param[in]     itypfl        indicator (take rho into account or not)
 *                               - 1 compute \f$ \rho\vect{u}\cdot\vect{s} \f$
 *                               - 0 compute \f$ \vect{u}\cdot\vect{s} \f$
 * \param[in]     iflmb0        the mass flux is set to 0 on walls and
 *                               symmetries if = 1
 * \param[in]     init          the mass flux is initialized to 0 if > 0
 * \param[in]     inc           indicator
 *                               - 0 solve an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     nswrgu        number of sweeps for the reconstruction
 *                               of the gradients
 * \param[in]     imligu        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thanks to neighbooring gradients
 *                               - = 1 thanks to the mean gradient
 * \param[in]     iwarnu        verbosity
 * \param[in]     epsrgu        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgu        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     rom           cell density
 * \param[in]     romb          density at boundary faces
 * \param[in]     vel           vector variable
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part - vector array )
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part - 3x3 tensor array)
 * \param[in,out] i_massflux    mass flux at interior faces \f$ \dot{m}_\fij \f$
 * \param[in,out] b_massflux    mass flux at boundary faces \f$ \dot{m}_\fib \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_mass_flux(const cs_mesh_t             *m,
             const cs_mesh_quantities_t  *fvq,
             int                          f_id,
             int                          itypfl,
             int                          iflmb0,
             int                          init,
             int                          inc,
             int                          imrgra,
             int                          nswrgu,
             int                          imligu,
             int                          iwarnu,
             double                       epsrgu,
             double                       climgu,
             const cs_real_t              rom[],
             const cs_real_t              romb[],
             const cs_real_3_t            vel[],
             const cs_real_3_t            coefav[],
             const cs_real_33_t           coefbv[],
             cs_real_t          *restrict i_massflux,
             cs_real_t          *restrict b_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integrated mass flux on the cells.
 *
 * \f[
 * \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_divergence(const cs_mesh_t          *m,
              int                       init,
              const cs_real_t           i_massflux[],
              const cs_real_t           b_massflux[],
              cs_real_t       *restrict diverg);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integrated mass flux on the cells for a tensor variable.
 *
 * \f[
 * \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_massflux    mass flux vector at interior faces
 * \param[in]     b_massflux    mass flux vector at boundary faces
 * \param[in,out] diverg        mass flux divergence vector
 */
/*----------------------------------------------------------------------------*/

void
cs_tensor_divergence(const cs_mesh_t            *m,
                     int                         init,
                     const cs_real_3_t           i_massflux[],
                     const cs_real_3_t           b_massflux[],
                     cs_real_3_t       *restrict diverg);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Project the external source terms to the faces in coherence with
 * cs_face_diffusion_scalar for the improved hydrostatic pressure algorithm
 * (iphydr=1).
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgu        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
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
cs_ext_force_flux(const cs_mesh_t          *m,
                  cs_mesh_quantities_t     *fvq,
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
/*!
 * \brief Project the external source terms to the faces in coherence with
 * cs_face_anisotropic_diffusion_scalar for the improved hydrostatic pressure
 * algorithm (iphydr=1).
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
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
cs_ext_force_anisotropic_flux(const cs_mesh_t          *m,
                              cs_mesh_quantities_t     *fvq,
                              int                       init,
                              int                       nswrgp,
                              int                       ircflp,
                              const cs_real_3_t         frcxt[],
                              const cs_real_t           cofbfp[],
                              const cs_real_t           i_visc[],
                              const cs_real_t           b_visc[],
                              cs_real_6_t               viscel[],
                              const cs_real_2_t         weighf[],
                              cs_real_t       *restrict i_massflux,
                              cs_real_t       *restrict b_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add \f$ \rho \tens{r} \vect{s}_\ij\f$ to a flux.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     f_id          field id (or -1)
 * \param[in]     itypfl        indicator (take rho into account or not)
 *                               - 1 compute \f$ \rho\vect{u}\cdot\vect{s} \f$
 *                               - 0 compute \f$ \vect{u}\cdot\vect{s} \f$
 * \param[in]     iflmb0        the mass flux is set to 0 on walls and
 *                               symmetries if = 1
 * \param[in]     init          the mass flux is initialized to 0 if > 0
 * \param[in]     inc           indicator
 *                               - 0 solve an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     nswrgu        number of sweeps for the reconstruction
 *                               of the gradients
 * \param[in]     imligu        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thanks to neighbooring gradients
 *                               - = 1 thanks to the mean gradient
 * \param[in]     iwarnu        verbosity
 * \param[in]     epsrgu        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgu        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     c_rho         cell density
 * \param[in]     b_rho         density at boundary faces
 * \param[in]     c_var         variable
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part - symmetric tensor array)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part - 6x6 symmetric tensor array)
 * \param[in,out] i_massflux    mass flux at interior faces \f$ \dot{m}_\fij \f$
 * \param[in,out] b_massflux    mass flux at boundary faces \f$ \dot{m}_\fib \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_tensor_face_flux(const cs_mesh_t             *m,
                    const cs_mesh_quantities_t  *fvq,
                    int                          f_id,
                    int                          itypfl,
                    int                          iflmb0,
                    int                          init,
                    int                          inc,
                    int                          imrgra,
                    int                          nswrgu,
                    int                          imligu,
                    int                          iwarnu,
                    double                       epsrgu,
                    double                       climgu,
                    const cs_real_t              c_rho[],
                    const cs_real_t              b_rho[],
                    const cs_real_6_t            c_var[],
                    const cs_real_6_t            coefav[],
                    const cs_real_66_t           coefbv[],
                    cs_real_3_t        *restrict i_massflux,
                    cs_real_3_t        *restrict b_massflux);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DIVERGENCE_H__ */
