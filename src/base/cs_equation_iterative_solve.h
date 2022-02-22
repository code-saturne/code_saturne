#ifndef __CS_EQUATION_ITERATIVE_SOLVE_H__
#define __CS_EQUATION_ITERATIVE_SOLVE_H__

/*============================================================================
 * Scalar convection diffusion equation solver.
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

*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_parameters.h"

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_equation_iterative_solve.c
 *
 * \brief This file gathers functions that solve advection diffusion equations
 * with source terms for one time step for a scalar, vector or tensor variable.
 *
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*! \brief This function solves an advection diffusion equation with source
 * terms for one time step for the variable \f$ a \f$.
 *
 * The equation reads:
 *
 * \f[
 * f_s^{imp}(a^{n+1}-a^n)
 * + \divs \left( a^{n+1} \rho \vect{u} - \mu \grad a^{n+1} \right)
 * = Rhs
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * f_s^{imp} \delta a
 * + \divs \left( \delta a \rho \vect{u} - \mu \grad \delta a \right)
 * = Rhs^1
 * \f]
 *
 * where \f$ \delta a = a^{n+1} - a^n\f$ and
 * \f$ Rhs^1 = Rhs - \divs( a^n \rho \vect{u} - \mu \grad a^n)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * f_s^{imp} \delta a^k
 * + \divs \left(\delta a^k \rho \vect{u}-\mu\grad\delta a^k \right)
 * = Rhs^k
 * \f]
 *
 * where \f$Rhs^k=Rhs-f_s^{imp}(a^k-a^n)
 * - \divs \left( a^k\rho\vect{u}-\mu\grad a^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ f_s^{imp} \f$ here!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     iterns        external sub-iteration number
 * \param[in]     f_id          field id (or -1)
 * \param[in]     iescap        compute the predictor indicator if 1
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     normp         Reference norm to solve the system (optional)
 *                              if negative: recomputed here
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ a^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ a^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
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
 * \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
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
 * \param[in]     rovsdt        \f$ f_s^{imp} \f$
 * \param[in]     smbrp         Right hand side \f$ Rhs^k \f$
 * \param[in,out] pvar          current variable
 * \param[out]    dpvar         last variable increment
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[out]    eswork        prediction-stage error estimator
 *                              (if iescap > 0)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_scalar(int                   idtvar,
                                   int                   iterns,
                                   int                   f_id,
                                   const char           *name,
                                   int                   iescap,
                                   int                   imucpp,
                                   cs_real_t             normp,
                                   cs_var_cal_opt_t     *var_cal_opt,
                                   const cs_real_t       pvara[],
                                   const cs_real_t       pvark[],
                                   const cs_real_t       coefap[],
                                   const cs_real_t       coefbp[],
                                   const cs_real_t       cofafp[],
                                   const cs_real_t       cofbfp[],
                                   const cs_real_t       i_massflux[],
                                   const cs_real_t       b_massflux[],
                                   const cs_real_t       i_viscm[],
                                   const cs_real_t       b_viscm[],
                                   const cs_real_t       i_visc[],
                                   const cs_real_t       b_visc[],
                                   cs_real_6_t           viscel[],
                                   const cs_real_2_t     weighf[],
                                   const cs_real_t       weighb[],
                                   int                   icvflb,
                                   const int             icvfli[],
                                   const cs_real_t       rovsdt[],
                                   cs_real_t             smbrp[],
                                   cs_real_t             pvar[],
                                   cs_real_t             dpvar[],
                                   const cs_real_t       xcpp[],
                                   cs_real_t             eswork[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function solves an advection diffusion equation with source terms
 * for one time step for the vector variable \f$ \vect{a} \f$.
 *
 * The equation reads:
 *
 * \f[
 * \tens{f_s}^{imp}(\vect{a}^{n+1}-\vect{a}^n)
 * + \divv \left( \vect{a}^{n+1} \otimes \rho \vect {u}
 *              - \mu \gradt \vect{a}^{n+1}\right)
 * = \vect{Rhs}
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \vect{a}
 * + \divv \left( \delta \vect{a} \otimes \rho \vect{u}
 *              - \mu \gradt \delta \vect{a} \right)
 * = \vect{Rhs}^1
 * \f]
 *
 * where \f$ \delta \vect{a} = \vect{a}^{n+1} - \vect{a}^n\f$ and
 * \f$ \vect{Rhs}^1 = \vect{Rhs}
 * - \divv \left( \vect{a}^n \otimes \rho \vect{u}
 *              - \mu \gradt \vect{a}^n \right)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \vect{a}^k
 * + \divv \left( \delta \vect{a}^k \otimes \rho \vect{u}
 *              - \mu \gradt \delta \vect{a}^k \right)
 * = \vect{Rhs}^k
 * \f]
 *
 * where \f$ \vect{Rhs}^k = \vect{Rhs}
 * - \tens{f_s}^{imp} \left(\vect{a}^k-\vect{a}^n \right)
 * - \divv \left( \vect{a}^k \otimes \rho \vect{u}
 *              - \mu \gradt \vect{a}^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     iterns        external sub-iteration number
 * \param[in]     f_id          field id (or -1)
 * \param[in]     name          associated name if f_id < 0, or NULL
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     iescap        compute the predictor indicator if 1
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafv        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfv        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     i_secvis      secondary viscosity at interior faces
 * \param[in]     b_secvis      secondary viscosity at boundary faces
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
 * \param[in,out] pvar          current variable
 * \param[out]    eswork        prediction-stage error estimator
 *                              (if iescap > 0)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_vector(int                   idtvar,
                                   int                   iterns,
                                   int                   f_id,
                                   const char           *name,
                                   int                   ivisep,
                                   int                   iescap,
                                   cs_var_cal_opt_t     *var_cal_opt,
                                   const cs_real_3_t     pvara[],
                                   const cs_real_3_t     pvark[],
                                   const cs_real_3_t     coefav[],
                                   const cs_real_33_t    coefbv[],
                                   const cs_real_3_t     cofafv[],
                                   const cs_real_33_t    cofbfv[],
                                   const cs_real_t       i_massflux[],
                                   const cs_real_t       b_massflux[],
                                   cs_real_t             i_viscm[],
                                   const cs_real_t       b_viscm[],
                                   const cs_real_t       i_visc[],
                                   const cs_real_t       b_visc[],
                                   const cs_real_t       secvif[],
                                   const cs_real_t       secvib[],
                                   cs_real_6_t           viscel[],
                                   const cs_real_2_t     weighf[],
                                   const cs_real_t       weighb[],
                                   int                   icvflb,
                                   const int             icvfli[],
                                   cs_real_33_t          fimp[],
                                   cs_real_3_t           smbrp[],
                                   cs_real_3_t           pvar[],
                                   cs_real_3_t           eswork[]);

/*----------------------------------------------------------------------------*/
/*! \brief This function solves an advection diffusion equation with source
 * terms for one time step for the symmetric tensor variable
 * \f$ \tens{\variat} \f$.
 *
 * The equation reads:
 *
 * \f[
 * \tens{f_s}^{imp}(\tens{\variat}^{n+1}-\tens{\variat}^n)
 * + \divt \left( \tens{\variat}^{n+1} \otimes \rho \vect {u}
 *              - \mu \gradtt \tens{\variat}^{n+1}\right)
 * = \tens{Rhs}
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \tens{\variat}
 * + \divt \left( \delta \tens{\variat} \otimes \rho \vect{u}
 *              - \mu \gradtt \delta \tens{\variat} \right)
 * = \tens{Rhs}^1
 * \f]
 *
 * where \f$ \delta \tens{\variat} = \tens{\variat}^{n+1} - \tens{\variat}^n\f$
 * and \f$ \tens{Rhs}^1 = \tens{Rhs}
 * - \divt \left( \tens{\variat}^n \otimes \rho \vect{u}
 *              - \mu \gradtt \tens{\variat}^n \right)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \tens{\variat}^k
 * + \divt \left( \delta \tens{\variat}^k \otimes \rho \vect{u}
 *              - \mu \gradtt \delta \tens{\variat}^k \right)
 * = \tens{Rhs}^k
 * \f]
 *
 * where \f$ \tens{Rhs}^k = \tens{Rhs}
 * - \tens{f_s}^{imp} \left(\tens{\variat}^k-\tens{\variat}^n \right)
 * - \divt \left( \tens{\variat}^k \otimes \rho \vect{u}
 *              - \mu \gradtt \tens{\variat}^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
 * \param[in]     coefats       boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbts       boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofafts       boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfts       boundary condition array for the diffusion
 *                               of the variable (Implicit part)
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
 * \param[in,out] pvar          current variable
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_tensor(int                   idtvar,
                                   int                   f_id,
                                   const char           *name,
                                   cs_var_cal_opt_t     *var_cal_opt,
                                   const cs_real_6_t     pvara[],
                                   const cs_real_6_t     pvark[],
                                   const cs_real_6_t     coefats[],
                                   const cs_real_66_t    coefbts[],
                                   const cs_real_6_t     cofafts[],
                                   const cs_real_66_t    cofbfts[],
                                   const cs_real_t       i_massflux[],
                                   const cs_real_t       b_massflux[],
                                   const cs_real_t       i_viscm[],
                                   const cs_real_t       b_viscm[],
                                   const cs_real_t       i_visc[],
                                   const cs_real_t       b_visc[],
                                   cs_real_6_t           viscel[],
                                   const cs_real_2_t     weighf[],
                                   const cs_real_t       weighb[],
                                   int                   icvflb,
                                   const int             icvfli[],
                                   const cs_real_66_t    fimp[],
                                   cs_real_6_t           smbrp[],
                                   cs_real_6_t           pvar[]);

END_C_DECLS

#endif /* __CS_EQUATION_ITERATIVE_SOLVE_H__ */
