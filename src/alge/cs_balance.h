#ifndef __CS_BALANCE_H__
#define __CS_BALANCE_H__

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion
 * terms of a transport equation of a scalar field \f$ \varia \f$.
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
 * - \f$ Rhs \f$ has already been initialized before calling bilsca!
 * - mind the minus sign
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: second order
 * - ischcp = 1: centered
 * - imucpp = 0: do not multiply the convective part by \f$ C_p \f$
 * - imucpp = 1: multiply the convective part by \f$ C_p \f$
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     iccocg        indicator
 *                               - 1 re-compute cocg matrix (for iterative gradients)
 *                               - 0 otherwise
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_scalar(int                idtvar,
                  int                f_id,
                  int                imucpp,
                  int                imasac,
                  int                inc,
                  int                iccocg,
                  cs_var_cal_opt_t  *var_cal_opt,
                  cs_real_t          pvar[],
                  const cs_real_t    pvara[],
                  const cs_real_t    coefap[],
                  const cs_real_t    coefbp[],
                  const cs_real_t    cofafp[],
                  const cs_real_t    cofbfp[],
                  const cs_real_t    i_massflux[],
                  const cs_real_t    b_massflux[],
                  const cs_real_t    i_visc[],
                  const cs_real_t    b_visc[],
                  cs_real_6_t        viscel[],
                  const cs_real_t    xcpp[],
                  const cs_real_2_t  weighf[],
                  const cs_real_t    weighb[],
                  int                icvflb,
                  const int          icvfli[],
                  cs_real_t          smbrp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion
 * terms of a transport equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilscv!
 * - mind the sign minus
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: second order
 * - ischcp = 1: centered
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in]     secvib        secondary viscosity at boundary faces
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbr          right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_vector(int                  idtvar,
                  int                  f_id,
                  int                  imasac,
                  int                  inc,
                  int                  ivisep,
                  cs_var_cal_opt_t    *var_cal_opt,
                  cs_real_3_t          pvar[],
                  const cs_real_3_t    pvara[],
                  const cs_real_3_t    coefav[],
                  const cs_real_33_t   coefbv[],
                  const cs_real_3_t    cofafv[],
                  const cs_real_33_t   cofbfv[],
                  const cs_real_t      i_massflux[],
                  const cs_real_t      b_massflux[],
                  const cs_real_t      i_visc[],
                  const cs_real_t      b_visc[],
                  const cs_real_t      secvif[],
                  const cs_real_t      secvib[],
                  cs_real_6_t          viscel[],
                  const cs_real_2_t    weighf[],
                  const cs_real_t      weighb[],
                  int                  icvflb,
                  const int            icvfli[],
                  cs_real_3_t          smbr[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion
 * terms of a transport equation of a tensor field \f$ \tens{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \tens{Rhs} = \tens{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \tens{\varia}_\fij - \tens{\varia}_\celli \right)
 *      - \mu_\fij \gradt_\fij \tens{\varia} \cdot \tens{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \tens{Rhs} \f$ has already been initialized before calling bilscts!
 * - mind the sign minus
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: second order
 * - ischcp = 1: centered
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     coefa       boundary condition array for the variable
 *                              (Explicit part)
 * \param[in]     coefb       boundary condition array for the variable
 *                              (Impplicit part)
 * \param[in]     cofaf       boundary condition array for the diffusion
 *                              of the variable (Explicit part)
 * \param[in]     cofbf       boundary condition array for the diffusion
 *                              of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_tensor(int                 idtvar,
                  int                 f_id,
                  int                 imasac,
                  int                 inc,
                  cs_var_cal_opt_t   *var_cal_opt,
                  cs_real_6_t         pvar[],
                  const cs_real_6_t   pvara[],
                  const cs_real_6_t   coefa[],
                  const cs_real_66_t  coefb[],
                  const cs_real_6_t   cofaf[],
                  const cs_real_66_t  cofbf[],
                  const cs_real_t     i_massflux[],
                  const cs_real_t     b_massflux[],
                  const cs_real_t     i_visc[],
                  const cs_real_t     b_visc[],
                  cs_real_6_t         viscel[],
                  const cs_real_2_t   weighf[],
                  const cs_real_t     weighb[],
                  int                 icvflb,
                  const int           icvfli[],
                  cs_real_6_t         smbrp[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BALANCE_H__ */
