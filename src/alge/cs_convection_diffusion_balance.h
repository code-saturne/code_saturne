#ifndef __CS_CONVECTION_DIFFUSION_BALANCE_H__
#define __CS_CONVECTION_DIFFUSION_BALANCE_H__

/*============================================================================
 * Explicit convection diffusion balance
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
 * Wrapper to cs_convection_diffusion_balance
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilsc5, BILSC5)
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
 cs_real_t                rhs[]);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*! \brief This function adds the explicit part of the convection/diffusion
  terms of a standard transport equation of a scalar field \f$ \varia \f$.

  More precisely, the right hand side \f$ Rhs \f$ is updated as
  follows:
  \f[
  Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
         \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
       - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
  \f]

  Warning:
  - \f$ Rhs \f$ has already been initialized before calling bilsc2!
  - mind the sign minus

  Options:
  - blencp = 0: upwind scheme for the advection
  - blencp = 1: no upwind scheme except in the slope test
  - ischcp = 0: second order
  - ischcp = 1: centred

*/
/*-------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 ______________________________________________________________________________*/
/*!
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
 *                               - 1 re-compute cocg matrix (for iterativ gradients)
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
 *                               (Impplicit part)
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
 * \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
 */
/*-------------------------------------------------------------------------------*/

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
                               cs_real_t *restrict       rhs);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CONVECTION_DIFFUSION_BALANCE_H__ */
