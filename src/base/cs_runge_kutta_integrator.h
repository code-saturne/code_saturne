#ifndef RK_INTEGRATOR_H
#define RK_INTEGRATOR_H

/*============================================================================
 * Explicit Runge-Kutta integrator utilities.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_balance.h"
#include "alge/cs_blas.h"

#include "base/cs_array.h"
#include "base/cs_base_accel.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "base/cs_dispatch.h"
#include "mesh/cs_mesh.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_runge_kutta_integrator_priv.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a RK integrator.
 *
 * \param[in] scheme  RK scheme type (RK_NONE, RK1, RK2, RK3, RK4)
 * \param[in] name    associated equation or field's name
 * \param[in] dt      time step
 * \param[in] dim     variable dimentsion
 * \param[in] n_elts  number of computational elements
 *
 * return the RK integrator's id in the RK list
 */
/*----------------------------------------------------------------------------*/

int
cs_runge_kutta_integrator_create(cs_runge_kutta_scheme_t   scheme,
                                 const char               *name,
                                 const cs_real_t          *dt,
                                 int                       dim,
                                 cs_lnum_t                 n_elts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create RK integrator structures
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrators_initialize();

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a RK integrator at the begining of a time step.
 *
 * \param[in]      ctx         Reference to dispatch context
 * \param[in,out]  rk          pointer to a Runge-Kutta integrator
 * \param[in]      rho         mass density within the geometry support
 * \param[in]      vol         measure of the geometry support
 * \param[in]      pvara       high level variable array at the begining
 *                             of a time step
 */
/*----------------------------------------------------------------------------*/

template<cs_lnum_t stride>
void
cs_runge_kutta_init_state(cs_dispatch_context          &ctx,
                          cs_runge_kutta_integrator_t  *rk,
                          const cs_real_t              *rho,
                          const cs_real_t              *vol,
                          cs_real_t                    *pvara)
{
  if (rk == nullptr)
    return;

  rk->i_stage = 0;

  const cs_lnum_t n_elts = rk->n_elts;

  cs_real_t *mass = rk->mass;
  cs_real_t *u_old = rk->u_old;

  ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t i_elt) {
    mass[i_elt] = rho[i_elt]*vol[i_elt];
    for (cs_lnum_t k = 0; k < stride; k++)
      u_old[stride*i_elt + k] = pvara[stride*i_elt + k];
  });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform one Runge-Kutta staging.
 *
 * \param[in]      ctx         Reference to dispatch context
 * \param[in,out]  rk          pointer to a Runge-Kutta integrator
 * \param[in,out]  pvar_stage  high level variable array along stages
 */
/*----------------------------------------------------------------------------*/

template<cs_lnum_t stride>
void
cs_runge_kutta_staging(cs_dispatch_context          &ctx,
                       cs_runge_kutta_integrator_t  *rk,
                       cs_real_t                    *pvar_stage)
{
  assert(rk != nullptr);

  const int n_elts = rk->n_elts;
  const cs_real_t *dt = rk->dt;
  const cs_real_t *mass = rk->mass;

  // get the current stage index
  const int i_stg = rk->i_stage;

  cs_real_t *u0    = rk->u_old;
  cs_real_t *u_new = rk->u_new;
  cs_real_t *rhs_stages = rk->rhs_stages;

  const cs_real_t *a = rk->rk_coeff.a + RK_HIGHEST_ORDER*i_stg;

  ctx.parallel_for (n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t i_elt) {
    for (cs_lnum_t k = 0; k < stride; k++)
      u_new[stride*i_elt + k] = u0[stride*i_elt + k];

    for (int j_stg = 0; j_stg <= i_stg; j_stg++) {
      cs_real_t *rhs  = rhs_stages + j_stg*stride*n_elts;
      for (cs_lnum_t k = 0; k < stride; k++)
        u_new[stride*i_elt + k] += dt[i_elt]/mass[i_elt]*a[j_stg]
                                 * rhs[i_elt*stride + k];
    }

    for (cs_lnum_t k = 0; k < stride; k++)
      pvar_stage[stride*i_elt + k] = u_new[stride*i_elt + k];
  });

  ctx.wait();

  rk->i_stage++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set initial rhs per stage for Runge-Kutta integrator.
 *        Align with legacy equations' building sequence.
 *        1) Collect initial rhs done in
 *        cs_solve_navier_stokes/cs_solve_equation_scalar
 *        2) Complete rhs with adding explicit part of the
 *        convection/diffusion balance
 *
 * \param[in]      ctx         Reference to dispatch context
 * \param[in,out]  rk          pointer to a Runge-Kutta integrator
 * \param[in]      rhs_pvar    pointer to high level rhs array
 */
/*----------------------------------------------------------------------------*/

template<int stride>
void
cs_runge_kutta_stage_set_initial_rhs(cs_dispatch_context         &ctx,
                                     cs_runge_kutta_integrator_t *rk,
                                     cs_real_t                   *rhs_pvar)
{
  assert(rk != nullptr);
  const int n_elts = rk->n_elts;
  // get the current stage index
  const int i_stg = rk->i_stage;

  cs_real_t *rhs  = rk->rhs_stages + i_stg*stride*n_elts;
  ctx.parallel_for (n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t i_elt) {
    for (cs_lnum_t k = 0; k < stride; k++)
      rhs[stride*i_elt + k] = rhs_pvar[stride*i_elt + k];
  });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief prepare and complete rhs per stage for Runge-Kutta integrator.
 *        Align with legacy equations' building sequence.
 *        1) Collect initial rhs done in
 *        cs_solve_navier_stokes/cs_solve_equation_scalar
 *        2) Complete rhs with adding explicit part of the
 *        convection/diffusion balance
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized
 *   before calling cs_balance_vector!
 * - mind the sign minus
 *
 * \param[in]      ctx           Reference to dispatch context
 * \param[in,out]  rk            pointer to a Runge-Kutta integrator
 * \param[in]      idtvar        indicator of the temporal scheme
 * \param[in]      f_id          field id (or -1)
 * \param[in]      imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]      eqp           pointer to a cs_equation_param_t structure which
 *                               contains variable calculation options
 * \param[in]      bc_coeffs     boundary condition structure for the variable
 * \param[in]      i_massflux    mass flux at interior faces
 * \param[in]      b_massflux    mass flux at boundary faces
 * \param[in]      i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                at interior faces for the r.h.s.
 * \param[in]      b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                                at boundary faces for the r.h.s.
 * \param[in]      viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]      weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]      weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]      icvflb        global indicator of boundary convection flux
 *                                - 0 upwind scheme at all boundary faces
 *                                - 1 imposed flux at some boundary faces
 * \param[in]      icvfli        boundary face indicator array of convection flux
 *                                - 0 upwind scheme
 *                                - 1 imposed flux
 * \param[in]      pvar          solved variable (current time step)
 * \param[in]      xcpp          array of specific heat (Cp)
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_stage_complete_scalar_rhs
  (cs_dispatch_context         &ctx,
   cs_runge_kutta_integrator_t *rk,
   int                          idtvar,
   int                          f_id,
   int                          imucpp,
   cs_equation_param_t         *eqp,
   cs_field_bc_coeffs_t        *bc_coeffs,
   const cs_real_t              i_massflux[],
   const cs_real_t              b_massflux[],
   const cs_real_t              i_visc[],
   const cs_real_t              b_visc[],
   cs_real_t                    viscel[][6],
   const cs_real_t              weighf[][2],
   const cs_real_t              weighb[],
   int                          icvflb,
   const int                    icvfli[],
   cs_real_t                    pvar[],
   const cs_real_t              xcpp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief prepare and complete rhs per stage for Runge-Kutta integrator.
 *        Align with legacy equations' building sequence.
 *        1) Collect initial rhs done in
 *        cs_solve_navier_stokes/cs_solve_equation_scalar
 *        2) Complete rhs with adding explicit part of the
 *        convection/diffusion balance
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
 * - \f$ \vect{Rhs} \f$ has already been initialized
 *   before calling cs_balance_vector!
 * - mind the sign minus
 *
 * \param[in]      ctx           Reference to dispatch context
 * \param[in,out]  rk            pointer to a Runge-Kutta integrator
 * \param[in]      idtvar        indicator of the temporal scheme
 * \param[in]      f_id          field id (or -1)
 * \param[in]      ivisep        indicator to take \f$ \divv
 *                                \left(\mu \gradt \transpose{\vect{a}} \right)
 *                                -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                                - 1 take into account,
 *                                - 0 otherwise
 * \param[in]      eqp           pointer to a cs_equation_param_t structure which
 *                               contains variable calculation options
 * \param[in]      bc_coeffs     boundary condition structure for the variable
 * \param[in]      i_massflux    mass flux at interior faces
 * \param[in]      b_massflux    mass flux at boundary faces
 * \param[in]      i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                at interior faces for the r.h.s.
 * \param[in]      b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                                at boundary faces for the r.h.s.
 * \param[in]      secvif        secondary viscosity at interior faces
 * \param[in]      secvib        secondary viscosity at boundary faces
 * \param[in]      viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]      weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]      weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]      icvflb        global indicator of boundary convection flux
 *                                - 0 upwind scheme at all boundary faces
 *                                - 1 imposed flux at some boundary faces
 * \param[in]      icvfli        boundary face indicator array of convection flux
 *                                - 0 upwind scheme
 *                                - 1 imposed flux
 * \param[in]      pvar          solved velocity (current time step)
 * \param[in]      rhs_pvar    pointer to high level rhs array
 */
/*----------------------------------------------------------------------------*/

template<int stride>
void
cs_runge_kutta_stage_complete_rhs(cs_dispatch_context         &ctx,
                                  cs_runge_kutta_integrator_t *rk,
                                  int                          idtvar,
                                  int                          f_id,
                                  int                          ivisep,
                                  cs_equation_param_t         *eqp,
                                  cs_field_bc_coeffs_t        *bc_coeffs,
                                  const cs_real_t              i_massflux[],
                                  const cs_real_t              b_massflux[],
                                  const cs_real_t              i_visc[],
                                  const cs_real_t              b_visc[],
                                  const cs_real_t              i_secvis[],
                                  const cs_real_t              b_secvis[],
                                  cs_real_t                    viscel[][6],
                                  const cs_real_t              weighf[][2],
                                  const cs_real_t              weighb[],
                                  int                          icvflb,
                                  const int                    icvfli[],
                                  cs_real_t                    pvar[][stride])
{
  assert(rk != nullptr);
  const int n_elts = rk->n_elts;
  // get the current stage index
  const int i_stg = rk->i_stage;

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t *rhs = rk->rhs_stages + i_stg*stride*n_elts;

  cs_alloc_mode_t amode = ctx.alloc_mode(false);

  /* Allocate non reconstructed face value only if presence of limiter */

  cs_field_t *f = nullptr;
  if (f_id > -1)
    f = cs_field_by_id(f_id);

  using var_t = cs_real_t[stride];

  /* We compute the total explicit balance. */

  int  inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
  int  imasac = 0;

  eqp->theta = 1;

  cs_boundary_conditions_update_bc_coeff_face_values_strided<stride>
    (ctx, f, bc_coeffs, inc, eqp, pvar);

  if (stride == 3)
    cs_balance_vector(idtvar,
                      f_id,
                      imasac,
                      inc,       /* inc */
                      ivisep,
                      eqp,
                      nullptr, /* pvar == pvara */
                      (const cs_real_3_t *)pvar,
                      bc_coeffs,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      i_secvis,
                      b_secvis,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      nullptr,
                      nullptr,
                      (cs_real_3_t *)rhs);

  else if (stride == 6)
    cs_balance_tensor(idtvar,
                      f_id,
                      imasac,
                      inc,       /* inc */
                      eqp,
                      nullptr, /* pvar == pvara */
                      (const cs_real_6_t *)pvar,
                      bc_coeffs,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      (cs_real_6_t *)rhs);

  eqp->theta = 0;

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if Runge-Kutta integrator is activated.
 * \param[in]  rk          pointer to a Runge-Kutta integrator
 */
/*----------------------------------------------------------------------------*/

inline bool
cs_runge_kutta_is_active(cs_runge_kutta_integrator_t  *rk)
{
  if (rk == nullptr)
    return false;
  else
    return (rk->scheme > CS_RK_NONE);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the scaling factor when conducting a projection per stage
 *
 * \param[in]  rk   pointer to a Runge-Kutta integrator
 */
/*----------------------------------------------------------------------------*/

inline cs_real_t *
cs_runge_kutta_get_projection_time_scale_by_stage
(
  cs_dispatch_context           &ctx,
  cs_runge_kutta_integrator_t  *rk
)
{
  assert(rk != nullptr);
  assert(rk->i_stage > 0);

  const cs_lnum_t n_elts = rk->n_elts;
  // get the current stage index
  const int i_stg = rk->i_stage;
  const cs_real_t *dt = rk->dt;

  const cs_real_t scaling_factor = rk->rk_coeff.c[i_stg - 1];
  cs_real_t *scaled_dt = rk->scaled_dt;

  ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
    scaled_dt[i] = dt[i] * scaling_factor;
  });
  ctx.wait();

  return rk->scaled_dt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief indicate if the Runge-Kutta integrator is staging.
 *
 * \param[in]  rk   pointer to a runge-kutta integrator
 */
/*----------------------------------------------------------------------------*/

inline bool
cs_runge_kutta_is_staging(cs_runge_kutta_integrator_t  *rk)
{
  if (rk == nullptr)
    return false;

  return (rk->i_stage < rk->n_stages);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a Runge-Kutta integrator by id.
 */
/*----------------------------------------------------------------------------*/

cs_runge_kutta_integrator_t *
cs_runge_kutta_integrator_by_id(int  rk_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clean all Runge-Kutta integrators.
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrators_destroy();

/*----------------------------------------------------------------------------*/

#endif /* __RK_INTEGRATOR_H__ */
