#ifndef __RK_INTEGRATOR_H__
#define __RK_INTEGRATOR_H__

/*============================================================================
 * Explicit Runge-Kutta integrator utilities.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "alge/cs_balance.h"
#include "alge/cs_blas.h"

#include "base/cs_array.h"
#include "base/cs_base_accel.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "base/cs_dispatch.h"
#include "base/cs_field_default.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_runge_kutta_integrator_priv.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Global variables
 *============================================================================*/

extern cs_runge_kutta_integrator_t **cs_glob_rk_lst = nullptr;
int    _n_rk_integrators = 0;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Fill integrator coeffs for a given scheme.
 * Reference: Accuracy analysis of explicit Runge–Kutta methods
 * applied to the incompressible Navier–Stokes equations.
 * J.Comput.Phys. 231 (2012) 3041–3063
 */
/*----------------------------------------------------------------------------*/

static void
_init_runge_kutta_scheme(cs_runge_kutta_integrator_t   *rk,
                         cs_runge_kutta_scheme_t        scheme)
{
  double *a = rk->rk_coeff.a;
  double *c = rk->rk_coeff.c;
  memset(a, 0., RK_HIGHEST_ORDER*RK_HIGHEST_ORDER*sizeof(double));
  memset(c, 0., RK_HIGHEST_ORDER*sizeof(double));

  switch (scheme) {

  case CS_RK_NONE:
    return;

  case CS_RK1:
    rk->n_stages = 1;
    a[0] = 1.;
    c[0] = 1.;
    break;

  case CS_RK2:
    rk->n_stages = 2;

    a[0] = 1.;
    a[4] = 0.5, a[5] = 0.5;

    c[0] = 1., c[1] = 1.;
    break;

  case CS_RK3:
    rk->n_stages = 3;

    a[0] = 1.0/3.0;
    a[4] = -1.0,    a[5] = 2.0;
    a[8] = 0.0,     a[9] = 3./4.,  a[10] = 1./4.;

    c[0] = 1.0/3.0, c[1] = 1.0, c[2] = 1.0;

    break;

  case CS_RK4:
    rk->n_stages = 4;

    a[0] = 1.0;
    a[4] = 3.0/8.0,  a[5] = 1.0/8.0;
    a[8] = -1.0/8.0, a[9] = -3.0/8.0,   a[10] = 3.0/2.0;
    a[12] = 1.0/6.0, a[13] = -1.0/18.0, a[14] = 2.0/3.0, a[15] = 2.0/9.0;
    c[0] = 1.0, c[1] = 1.0/2.0, c[2] = 1.0, c[3] = 1.0;

    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
            "%s: Type of Runge-Kutta not available.\n"
            "%s: Stop building Runge-Kutta integrator.\n",
            __func__, __func__);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a RK integrator.
 *
 * \param[in] scheme  RK scheme type (RK_NONE, RK1, RK2, RK3, RK4)
 * \param[in] name    associated equation or field's name
 * \param[in] dt      time step
 * \param[in] n_elts  number of computational elements
 *
 * return the RK integrator's id in the RK list
 */
/*----------------------------------------------------------------------------*/
template<int stride>
static int
_runge_kutta_create(cs_runge_kutta_scheme_t      scheme,
                    const char                  *name,
                    const cs_real_t             *dt,
                    const cs_lnum_t              n_elts)
{
  cs_runge_kutta_integrator_t *rk = nullptr;

  CS_MALLOC(rk, 1, cs_runge_kutta_integrator_t);

  rk->scheme = scheme;

  CS_MALLOC(rk->name, strlen(name) + 1, char);
  strcpy(rk->name, name);

  rk->dt = dt;

  rk->n_elts = n_elts;
  rk->n_stages = 1;
  rk->i_stage = 0;

  CS_MALLOC_HD(rk->rk_coeff.a, RK_HIGHEST_ORDER*RK_HIGHEST_ORDER,
      double, cs_alloc_mode);
  CS_MALLOC_HD(rk->rk_coeff.c, RK_HIGHEST_ORDER,
      double, cs_alloc_mode);

  _init_runge_kutta_scheme(rk, scheme);

  CS_MALLOC_HD(rk->rhs_stages, rk->n_stages*stride*n_elts,
      cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(rk->u_old, stride*n_elts, cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(rk->scaled_dt, n_elts, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(rk->mass, n_elts, cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(rk->u_new, stride*n_elts, cs_real_t, cs_alloc_mode);

  int _rk_integrator_id = _n_rk_integrators;

  CS_REALLOC(cs_glob_rk_lst, _n_rk_integrators + 1,
            cs_runge_kutta_integrator_t*);

  cs_glob_rk_lst[_rk_integrator_id] = rk;

  _n_rk_integrators++;

  return _rk_integrator_id;
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create RK integrator structures
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrators_initialize()
{
  const int n_fields = cs_field_n_fields();
  const cs_real_t *dt = CS_F_(dt)->val;

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      if (eqp->rk_def.scheme <= CS_RK_NONE)
        continue;

      switch (eqp->dim) {
        case 1:
          eqp->rk_def.rk_id =
            _runge_kutta_create<1>(eqp->rk_def.scheme, f->name, dt,
                n_cells_ext);
          break;

        case 3:
          eqp->rk_def.rk_id =
            _runge_kutta_create<3>(eqp->rk_def.scheme, f->name, dt,
                n_cells_ext);
          break;

        case 6:
          eqp->rk_def.rk_id =
            _runge_kutta_create<6>(eqp->rk_def.scheme, f->name, dt,
                n_cells_ext);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
              "%s: Equation's dimension for Runge-Kutta not available.\n"
              "%s: Stop building Runge-Kutta integrator.\n",
              __func__, __func__);
          break;
      }
    }
  }
}

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
template<int stride>
void
cs_runge_kutta_init_state(cs_dispatch_context         &ctx,
                          cs_runge_kutta_integrator_t *rk,
                          const cs_real_t             *rho,
                          const cs_real_t             *vol,
                          cs_real_t                   *pvara)
{
  if (rk == nullptr)
    return;

  rk->i_stage = 0;

  const cs_lnum_t n_elts = rk->n_elts;

  cs_real_t *mass = rk->mass;
  cs_real_t *u_old = rk->u_old;

  ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t i_elt) {
    mass[i_elt] = rho[i_elt]*vol[i_elt];
    for (int k = 0; k < stride; k++)
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
template<int stride>
void
cs_runge_kutta_staging(cs_dispatch_context         &ctx,
                       cs_runge_kutta_integrator_t *rk,
                       cs_real_t                   *pvar_stage)
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
    for (int k = 0; k < stride; k++)
      u_new[stride*i_elt + k] = u0[stride*i_elt + k];

    for (int j_stg = 0; j_stg <= i_stg; j_stg++) {
      cs_real_t *rhs  = rhs_stages + j_stg*stride*n_elts;
      for (int k = 0; k < stride; k++)
        u_new[stride*i_elt + k] += dt[i_elt]/mass[i_elt]*a[j_stg]
                                 * rhs[i_elt*stride + k];
    }

    for (int k = 0; k < stride; k++)
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
    for (int k = 0; k < stride; k++)
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

  int df_limiter_id = -1;
  cs_field_t *f = nullptr;
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
  }

  const int ircflp = eqp->ircflu;
  const int ircflb = (ircflp > 0) ? eqp->b_diff_flux_rc : 0;

  cs_bc_coeffs_solve_t bc_coeffs_solve;
  cs_init_bc_coeffs_solve(bc_coeffs_solve,
                          n_b_faces,
                          stride,
                          amode,
                          (df_limiter_id > -1 || ircflb != 1));

  using var_t = cs_real_t[stride];

  var_t *val_ip = (var_t *)bc_coeffs_solve.val_ip;
  var_t *val_f = (var_t *)bc_coeffs_solve.val_f;
  var_t *flux =  (var_t *)bc_coeffs_solve.flux;
  var_t *flux_lim =  (var_t *)bc_coeffs_solve.flux_lim;

  /* We compute the total explicit balance. */

  int  inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
  int  imasac = 0;

  eqp->theta = 1;

  cs_boundary_conditions_update_bc_coeff_face_values_strided<stride>
    (ctx, f, bc_coeffs, inc, eqp, pvar,
     val_ip, val_f, flux, flux_lim);

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
                      (const cs_real_3_t *)val_f,
                      (const cs_real_3_t *)flux_lim,
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
                      (const cs_real_6_t *)val_f,
                      (const cs_real_6_t *)flux_lim,
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

  cs_clear_bc_coeffs_solve(bc_coeffs_solve);
  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if Runge-Kutta integrator is activated.
 * \param[in]  rk          pointer to a Runge-Kutta integrator
 */
/*----------------------------------------------------------------------------*/

bool
cs_runge_kutta_is_active(cs_runge_kutta_integrator_t *rk)
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
 * \param[in]  rk          pointer to a Runge-Kutta integrator
 * */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_runge_kutta_get_projection_time_scale_by_stage
(cs_runge_kutta_integrator_t *rk)
{
  assert(rk != nullptr);
  assert(rk->i_stage > 0);

  const cs_lnum_t n_elts = rk->n_elts;
  // get the current stage index
  const int i_stg = rk->i_stage;
  const cs_real_t *dt = rk->dt;

  const cs_real_t scaling_factor = rk->rk_coeff.c[i_stg - 1];

  memset(rk->scaled_dt, 0., n_elts*sizeof(cs_real_t));

  cs_axpy(n_elts, scaling_factor, dt, rk->scaled_dt);

  return rk->scaled_dt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief indicate if the Runge-Kutta integrator is staging
 *
 * \param[in]  rk          pointer to a runge-kutta integrator
 * */
/*----------------------------------------------------------------------------*/

bool
cs_runge_kutta_is_staging(cs_runge_kutta_integrator_t *rk)
{
  if (rk == nullptr)
    return false;

  return (rk->i_stage < rk->n_stages);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a Runge-Kutta integrator by id
 * */
/*----------------------------------------------------------------------------*/

cs_runge_kutta_integrator_t *
cs_runge_kutta_integrator_by_id(int     rk_id)
{
  if (rk_id < 0
   || rk_id >= _n_rk_integrators)
    return nullptr;

  return cs_glob_rk_lst[rk_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory
 * \param[in]  rk      double pointer to a runge-kutta integrator
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrator_free(cs_runge_kutta_integrator_t **rk)
{
  cs_runge_kutta_integrator_t *_rk = *rk;

  if (_rk == nullptr)
    return;

  CS_FREE_HD(_rk->rhs_stages);
  CS_FREE_HD(_rk->rk_coeff.a);
  CS_FREE_HD(_rk->rk_coeff.c);

  CS_FREE(_rk->name);

  CS_FREE_HD(_rk->scaled_dt);
  CS_FREE_HD(_rk->mass);

  CS_FREE_HD(_rk->u_old);
  CS_FREE_HD(_rk->u_new);
  CS_FREE(_rk);

  *rk = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clean all Runge-Kutta integrators
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrators_destroy()
{
  for (int i = 0; i < _n_rk_integrators; i++)
    cs_runge_kutta_integrator_free(&cs_glob_rk_lst[i]);

  CS_FREE(cs_glob_rk_lst);
  cs_glob_rk_lst = nullptr;
  _n_rk_integrators = 0;
}

/*----------------------------------------------------------------------------*/

#endif /* __RK_INTEGRATOR_H__ */
