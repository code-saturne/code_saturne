#ifndef __CS_CONVECTION_DIFFUSION_H__
#define __CS_CONVECTION_DIFFUSION_H__

/*============================================================================
 * Convection-diffusion operators.
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
 * Wrapper to cs_convection_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc2, BILSC2)
(
 const cs_int_t  *const   idtvar,
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   ischcp,
 const cs_int_t  *const   isstpp,
 const cs_int_t  *const   icvflb,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   ifaccp,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   blencp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 const cs_real_t *const   relaxp,
 const cs_real_t *const   thetap,
 cs_real_t                pvar[],
 const cs_real_t          pvara[],
 const cs_int_t           bc_type[],
 const cs_int_t           icvfli[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 cs_real_t                rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc4, BILSC4)
(
 const cs_int_t  *const   idtvar,
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   ischcp,
 const cs_int_t  *const   isstpp,
 const cs_int_t  *const   icvflb,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   ifaccp,
 const cs_int_t  *const   ivisep,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   blencp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   relaxp,
 const cs_real_t *const   thetap,
 cs_real_3_t              pvar[],
 const cs_real_3_t        pvara[],
 const cs_int_t           bc_type[],
 const cs_int_t           icvfli[],
 const cs_real_3_t        coefav[],
 const cs_real_33_t       coefbv[],
 const cs_real_3_t        cofafv[],
 const cs_real_33_t       cofbfv[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          secvif[],
 cs_real_3_t              rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_convection_diffusion_thermal
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsct, BILSCT)
(
 const cs_int_t  *const   idtvar,
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   iconvp,
 const cs_int_t  *const   idiffp,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   ischcp,
 const cs_int_t  *const   isstpp,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   ifaccp,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   blencp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   extrap,
 const cs_real_t *const   relaxp,
 const cs_real_t *const   thetap,
 cs_real_t                pvar[],
 const cs_real_t          pvara[],
 const cs_int_t           bc_type[],
 const cs_real_t          coefap[],
 const cs_real_t          coefbp[],
 const cs_real_t          cofafp[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          xcpp[],
 cs_real_t                rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (diftnv, DIFTNV)
(
 const cs_int_t  *const   idtvar,
 const cs_int_t  *const   f_id,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   ifaccp,
 const cs_int_t  *const   ivisep,
 const cs_int_t  *const   iwarnp,
 const cs_real_t *const   epsrgp,
 const cs_real_t *const   climgp,
 const cs_real_t *const   relaxp,
 const cs_real_t *const   thetap,
 cs_real_3_t              pvar[],
 const cs_real_3_t        pvara[],
 const cs_int_t           bc_type[],
 const cs_real_3_t        coefav[],
 const cs_real_33_t       coefbv[],
 const cs_real_3_t        cofafv[],
 const cs_real_33_t       cofbfv[],
 const cs_real_33_t       i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          secvif[],
 cs_real_3_t              rhs[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_diffusion_scalar
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
 cs_real_t                b_massflux[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_anisotropic_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrmav, ITRMAV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
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
 cs_real_6_t              viscel[],
 const cs_real_2_t        weighf[],
 const cs_real_t          weighb[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_diffusion_scalar
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
 cs_real_t                diverg[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_anisotropic_diffusion_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (itrgrv, ITRGRV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   inc,
 const cs_int_t  *const   imrgra,
 const cs_int_t  *const   iccocg,
 const cs_int_t  *const   nswrgp,
 const cs_int_t  *const   imligp,
 const cs_int_t  *const   ircflp,
 const cs_int_t  *const   iphydp,
 const cs_int_t  *const   iwarnp,
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
 cs_real_6_t              viscel[],
 const cs_real_2_t        weighf[],
 const cs_real_t          weighb[],
 cs_real_t                diverg[]
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before calling bilsc2!
 * - mind the sign minus
 *
 * Options:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: second order
 * - ischcp = 1: centred
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     iconvp        indicator
 *                               - 1 convection,
 *                               - 0 sinon
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion,
 *                               - 0 sinon
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     ischcp        indicator
 *                               - 1 centred
 *                               - 0 2nd order
 * \param[in]     isstpp        indicator
 *                               - 1 without slope test
 *                               - 0 with slope test
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix
 *                                   (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     iwarnp        verbosity
 * \param[in]     blencp        fraction of upwinding
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     relaxp        coefficient of relaxation
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_scalar(
                               int                       idtvar,
                               int                       f_id,
                               int                       iconvp,
                               int                       idiffp,
                               int                       nswrgp,
                               int                       imligp,
                               int                       ircflp,
                               int                       ischcp,
                               int                       isstpp,
                               int                       icvflb,
                               int                       inc,
                               int                       imrgra,
                               int                       iccocg,
                               int                       ifaccp,
                               int                       iwarnp,
                               double                    blencp,
                               double                    epsrgp,
                               double                    climgp,
                               double                    extrap,
                               double                    relaxp,
                               double                    thetap,
                               cs_real_t       *restrict pvar,
                               const cs_real_t *restrict pvara,
                               const cs_int_t            bc_type[],
                               const cs_int_t            icvfli[],
                               const cs_real_t           coefap[],
                               const cs_real_t           coefbp[],
                               const cs_real_t           cofafp[],
                               const cs_real_t           cofbfp[],
                               const cs_real_t           i_massflux[],
                               const cs_real_t           b_massflux[],
                               const cs_real_t           i_visc[],
                               const cs_real_t           b_visc[],
                               cs_real_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * Options:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: second order
 * - ischcp = 1: centred
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     iconvp        indicator
 *                               - 1 convection,
 *                               - 0 sinon
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion,
 *                               - 0 sinon
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     ischcp        indicator
 *                               - 1 centred
 *                               - 0 2nd order
 * \param[in]     isstpp        indicator
 *                               - 1 without slope test
 *                               - 0 with slope test
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     iwarnp        verbosity
 * \param[in]     blencp        fraction of upwinding
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     relaxp        coefficient of relaxation
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     coefav        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_vector(
                               int                         idtvar,
                               int                         f_id,
                               int                         iconvp,
                               int                         idiffp,
                               int                         nswrgp,
                               int                         imligp,
                               int                         ircflp,
                               int                         ischcp,
                               int                         isstpp,
                               int                         icvflb,
                               int                         inc,
                               int                         imrgra,
                               int                         ifaccp,
                               int                         ivisep,
                               int                         iwarnp,
                               double                      blencp,
                               double                      epsrgp,
                               double                      climgp,
                               double                      relaxp,
                               double                      thetap,
                               cs_real_3_t       *restrict pvar,
                               const cs_real_3_t *restrict pvara,
                               const cs_int_t              bc_type[],
                               const cs_int_t              icvfli[],
                               const cs_real_3_t           coefav[],
                               const cs_real_33_t          coefbv[],
                               const cs_real_3_t           cofafv[],
                               const cs_real_33_t          cofbfv[],
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               const cs_real_t             secvif[],
                               cs_real_3_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 * equation of a scalar field \f$ \varia \f$ such as the temperature.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs + \sum_{\fij \in \Facei{\celli}}      \left(
 *        C_p\dot{m}_\ij \varia_\fij
 *      - \lambda_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * \f$ Rhs \f$ has already been initialized before calling bilsct!
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: second order
 * - ischcp = 1: centred
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     iconvp        indicator
 *                               - 1 convection,
 *                               - 0 sinon
 * \param[in]     idiffp        indicator
 *                               - 1 diffusion,
 *                               - 0 sinon
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     ischcp        indicator
 *                               - 1 centred
 *                               - 0 2nd order
 * \param[in]     isstpp        indicator
 *                               - 1 without slope test
 *                               - 0 with slope test
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
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     iwarnp        verbosity
 * \param[in]     blencp        fraction of upwinding
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     relaxp        coefficient of relaxation
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     coefap        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$)
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_thermal(
                                int                       idtvar,
                                int                       f_id,
                                int                       iconvp,
                                int                       idiffp,
                                int                       nswrgp,
                                int                       imligp,
                                int                       ircflp,
                                int                       ischcp,
                                int                       isstpp,
                                int                       inc,
                                int                       imrgra,
                                int                       iccocg,
                                int                       ifaccp,
                                int                       iwarnp,
                                double                    blencp,
                                double                    epsrgp,
                                double                    climgp,
                                double                    extrap,
                                double                    relaxp,
                                double                    thetap,
                                cs_real_t       *restrict pvar,
                                const cs_real_t *restrict pvara,
                                const cs_int_t            bc_type[],
                                const cs_real_t           coefap[],
                                const cs_real_t           coefbp[],
                                const cs_real_t           cofafp[],
                                const cs_real_t           cofbfp[],
                                const cs_real_t           i_massflux[],
                                const cs_real_t           b_massflux[],
                                const cs_real_t           i_visc[],
                                const cs_real_t           b_visc[],
                                const cs_real_t           xcpp[],
                                cs_real_t       *restrict rhs);

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensorial
 * diffusivity for a transport equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling diftnv!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     ifaccp        indicator
 *                               - 1 coupling activated
 *                               - 0 coupling not activated
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     relaxp        coefficient of relaxation
 * \param[in]     thetap        weightening coefficient for the theta-schema,
 *                               - thetap = 0: explicit scheme
 *                               - thetap = 0.5: time-centred
 *                               scheme (mix between Crank-Nicolson and
 *                               Adams-Bashforth)
 *                               - thetap = 1: implicit scheme
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_type       boundary condition type
 * \param[in]     coefav        boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (Implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_visc        \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf\jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_vector(
                                int                         idtvar,
                                int                         f_id,
                                int                         nswrgp,
                                int                         imligp,
                                int                         ircflp,
                                int                         inc,
                                int                         imrgra,
                                int                         ifaccp,
                                int                         ivisep,
                                int                         iwarnp,
                                double                      epsrgp,
                                double                      climgp,
                                double                      relaxp,
                                double                      thetap,
                                cs_real_3_t       *restrict pvar,
                                const cs_real_3_t *restrict pvara,
                                const cs_int_t              bc_type[],
                                const cs_real_3_t           coefav[],
                                const cs_real_33_t          coefbv[],
                                const cs_real_3_t           cofafv[],
                                const cs_real_33_t          cofbfv[],
                                const cs_real_33_t          i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t             secvif[],
                                cs_real_3_t       *restrict rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the face mass flux with the face pressure (or pressure
 * increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
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
 *                               (Implicit part)
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
cs_face_diffusion_scalar(const cs_mesh_t          *m,
                         cs_mesh_quantities_t     *fvq,
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
/*!
 * \brief Add the explicit part of the pressure gradient term to the mass flux
 * in case of anisotropic diffusion of the pressure field \f$ P \f$.
 *
 * More precisely, the mass flux side \f$ \dot{m}_\fij \f$ is updated as
 * follows:
 * \f[
 * \dot{m}_\fij = \dot{m}_\fij -
 *              \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
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
                                    (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
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
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_diffusion_scalar(const cs_mesh_t          *m,
                                     cs_mesh_quantities_t     *fvq,
                                     int                       init,
                                     int                       inc,
                                     int                       imrgra,
                                     int                       iccocg,
                                     int                       nswrgp,
                                     int                       imligp,
                                     int                       ircflp,
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
                                     cs_real_6_t     *restrict viscel,
                                     const cs_real_2_t         weighf[],
                                     const cs_real_t           weighb[],
                                     cs_real_t       *restrict i_massflux,
                                     cs_real_t       *restrict b_massflux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the cell mass flux divergence with the face pressure (or
 * pressure increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \sum_j \Delta t \grad_\fij p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
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
 *                               (Implicit part)
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
cs_diffusion_scalar(const cs_mesh_t          *m,
                    cs_mesh_quantities_t     *fvq,
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
/*!
 * \brief Add the explicit part of the divergence of the mass flux due to the
 * pressure gradient (routine analog to diften).
 *
 * More precisely, the divergence of the mass flux side
 * \f$ \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij \f$ is updated as follows:
 * \f[
 * \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  = \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  - \sum_{\fij \in \Facei{\celli}}
 *    \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
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
                                     (for iterativ gradients)
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     extrap        coefficient for extrapolation of the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
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
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] diverg        divergence of the mass flux
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_scalar(const cs_mesh_t          *m,
                                cs_mesh_quantities_t     *fvq,
                                int                       init,
                                int                       inc,
                                int                       imrgra,
                                int                       iccocg,
                                int                       nswrgp,
                                int                       imligp,
                                int                       ircflp,
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
                                cs_real_6_t     *restrict viscel,
                                const cs_real_2_t         weighf[],
                                const cs_real_t           weighb[],
                                cs_real_t       *restrict diverg);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CONVECTION_DIFFUSION_H__ */
