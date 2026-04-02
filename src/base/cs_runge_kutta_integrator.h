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

#if defined(__cplusplus)
namespace cs {

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic Runge-Kutta integrator class
 */
/*----------------------------------------------------------------------------*/

class runge_kutta_integrator {

/*=============================================================================
 * Public methods
 *============================================================================*/

public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  runge_kutta_integrator()
  {}

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor method with input arguments.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  runge_kutta_integrator
  (
    cs_runge_kutta_scheme_t  scheme, /*!<[in] Selected RK scheme */
    const char              *name,   /*!<[in] Associated equation or field's
                                              name */
    const cs_real_t         *dt,     /*!<[in] time step array (ptr) */
    int                      dim,    /*!<[in] Field's dimension */
    cs_lnum_t                n_elts  /*!<[in] Number of computational elements */
  )
  {
    _scheme = scheme;
    CS_MALLOC(_name, strlen(name) + 1, char);
    strcpy(_name, name);

    _dt.update_data((cs_real_t *)dt, n_elts);

    _stride = dim;
    _n_elts = n_elts;
    _n_stages = 1;
    _i_stage = 0;

    _a.reshape(RK_HIGHEST_ORDER, RK_HIGHEST_ORDER);
    _c.reshape(RK_HIGHEST_ORDER);

    init_scheme_();

    _rhs_stages.reshape(_n_stages, _n_elts, _stride);
    _u_old.reshape(_n_elts, _stride);
    _u_new.reshape(_n_elts, _stride);
    _mass.reshape(_n_elts);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  ~runge_kutta_integrator()
  {
    CS_FREE(_name);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief IDefault constructor method.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  void
  init_state
  (
    cs_dispatch_context &ctx,  /*!<[in] Reference to dispatch context */
    const cs_real_t     *rho,  /*!<[in] mass density within the geometry
                                        support */
    const cs_real_t     *vol,  /*!<[in] measure of the geometry support */
    cs_real_t           *pvara /*!<[in] high level variable array at the
                                        begining of a time step */
  )
  {
    _i_stage = 0;

    auto mass = _mass.view();
    auto u_old = _u_old.view();

    ctx.parallel_for(_n_elts, CS_LAMBDA (cs_lnum_t e_id) {
      mass[e_id] = rho[e_id] * vol[e_id];
      for (cs_lnum_t i = 0; i < _stride; i++)
        u_old(e_id, i) = pvara[_stride*e_id + i];
    });
    ctx.wait();
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Perform one Runge-Kutta staging.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  void
  solve_stage
  (
    cs_dispatch_context &ctx,       /*!<[in] Reference to dispatch context */
    cs_real_t           *pvar_stage /*!<[in] high level variable array along
                                             stages */
  )
  {
    auto dt = _dt.view_1d();
    auto mass = _mass.view();

    const int i_stg = _i_stage;

    auto u0 = _u_old.view();
    auto u_new = _u_new.view();
    auto rhs_stage = _rhs_stages.view();

    const auto a = _a.sub_view(i_stg);

    ctx.parallel_for(_n_elts, CS_LAMBDA (cs_lnum_t e_id) {
      for (cs_lnum_t i = 0; i < _stride; i++)
        u_new(e_id, i) = u0(e_id, i);

      for (int j_stg = 0; j_stg <= i_stg; j_stg++) {
        for (cs_lnum_t i = 0; i < _stride; i++)
          u_new(e_id, i) += rhs_stage(j_stg, e_id, i)* a[j_stg] * dt[e_id]
                          / mass[e_id];
      }

      for (cs_lnum_t i = 0; i < _stride; i++)
        pvar_stage[_stride*e_id + i] = u_new(e_id, i);
    });

    ctx.wait();

    _i_stage += 1;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Perform one Runge-Kutta staging for the potential in the NS system.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  void
  stage_projection_rhs
  (
    cs_dispatch_context &ctx,    /*!<[in] Reference to dispatch context */
    cs_real_t           *grad_dp /*!<[in] pressure increment gradient */
  )
  {
    const int i_stg = _i_stage - 1;
    const auto a = _a.sub_view(i_stg);

    auto rhs = _rhs_stages.sub_view(i_stg);
    auto mass = _mass.view();

    ctx.parallel_for(_n_elts, CS_LAMBDA (cs_lnum_t e_id) {
      for (cs_lnum_t i = 0; i < _stride; i++) {
        rhs(e_id, i) -= grad_dp[_stride*e_id + i] * mass[e_id];
        grad_dp[_stride*e_id + i] *= a[i_stg];
      }
    });
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set initial rhs per stage for Runge-Kutta integrator.
   *        Align with legacy equations' building sequence.
   *        1) Collect initial rhs done in
   *        cs_solve_navier_stokes/cs_solve_equation_scalar
   *        2) Complete rhs with adding explicit part of the
   *        convection/diffusion balance
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  void
  stage_set_initial_rhs
  (
    cs_dispatch_context &ctx,     /*!<[in] Reference to dispatch context */
    cs_real_t           *rhs_pvar /*!<[in] pointer to high level rhs array */
  )
  {
    const int i_stg = _i_stage - 1;
    auto rhs = _rhs_stages.sub_view(i_stg);

    ctx.parallel_for(_n_elts, CS_LAMBDA (cs_lnum_t e_id) {
      for (cs_lnum_t i = 0; i < _stride; i++)
        rhs(e_id, i) = rhs_pvar[_stride*e_id + i];
    });
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Check if integrator is active.
   *
   * \return true if active, false otherwise.
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline bool
  is_active() const
  {
    return (_scheme > CS_RK_NONE);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Check if integrator is still in the staging process
   *
   * \return true if staging, false otherwise
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline bool
  is_staging() const
  {
    return (_i_stage < _n_stages);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get scheme used by integrator
   *
   * \return sheme type of the integrator (cs_runge_kutta_scheme_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline cs_runge_kutta_scheme_t
  scheme() const
  {
    return _scheme;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get number of elements used for resolution.
   *
   * \return number of elements (cs_lnum_t)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline cs_lnum_t
  n_elts() const
  {
    return _n_elts;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get current stage of RK resolution
   *
   * \return current stage id (int)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline int
  i_stage() const
  {
    return _i_stage;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get solved variable stride
   *
   * \return variable stride (int)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST_DEVICE
  inline int
  stride() const
  {
    return _stride;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get raw pointer to rhs of a given stage
   *
   * \return raw pointer (cs_real_t *)
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  cs_real_t *
  get_rhs_stage_sub_array
  (
    int stage_id /*!<[in] stage id */
  )
  {
    return _rhs_stages.sub_array(stage_id);
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get "a" coefficients for a given stage
   *
   * \return 1D span (cs_span<cs_real_t>) of the coefficients
   */
  /*--------------------------------------------------------------------------*/
  CS_F_HOST
  cs_span<cs_real_t>
  get_stage_coeff_a
  (
    int i /*!<[in] stage id */
  )
  {
    return _a.sub_view(i);
  }

private:

  /*===========================================================================
   * Private methods
   *==========================================================================*/

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief  Fill integrator coeffs for a given scheme.
   *
   * Reference: Accuracy analysis of explicit Runge–Kutta methods
   * applied to the incompressible Navier–Stokes equations.
   * J.Comput.Phys. 231 (2012) 3041–3063
   */
  /*--------------------------------------------------------------------------*/

  CS_F_HOST
  void
  init_scheme_()
  {
    _a.zero();
    _c.zero();

    switch (_scheme) {
    case CS_RK_NONE:
      return;

    case CS_RK1:
      _n_stages = 1;
      _a(0, 0) = 1.;
      _c[0]    = 1.;

      break;

    case CS_RK2:
      [[fallthrough]];
    case CS_RK2_HEUN:
      _n_stages = 2;
      _a(0,0) = 1.;
      _a(1,0) = 0.5, _a(1,1) = 0.5;

      _c[0] = 1.;
      _c[1] = 1.;

      break;

    case CS_RK3:
      [[fallthrough]];
    case CS_RK3_WRAY:
      _n_stages = 3;

      _a(0,0) = 8. / 15.;
      _a(1,0) = 1./4., _a(1,1) = 5./12.;
      _a(2,0) = 1./4., _a(2,1) = 0., _a(2,2) = 3./4.;

      _c[0] = 8./15.;
      _c[1] = 2./3.;
      _c[2] = 1.;

      break;

    case CS_RK3_SSP:
      _n_stages = 3;
      _a(0,0) = 1.;
      _a(1,0) = 1./4., _a(1,1) = 1./4.;
      _a(2,0) = 1./6., _a(2,1) = 1./6., _a(2,2) = 2./3.;

      _c[0] = 1.;
      _c[1] = 1./2.;
      _c[2] = 1.;

      break;

    case CS_RK4:
      _n_stages = 4;
      _a(0,0) = 1.;
      _a(1,0) = 3./8., _a(1,1) = 1./8.;
      _a(2,0) = -1./8., _a(2,1) = -3./8., _a(2,2) = 3./2.;
      _a(3,0) = 1./6., _a(3,1) = -1./18., _a(3,2) = 2./3., _a(3,3) = 2./9.;

      _c[0] = 1.;
      _c[1] = 1./2.;
      _c[2] = 1.;
      _c[3] = 1.;

      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Type of Runge-Kutta not available.\n"
                "%s: Stop building Runge-Kutta integrator.\n",
                __func__, __func__);
      break;
    }
  }

  /*===========================================================================
   * Private members
   *==========================================================================*/

  cs_runge_kutta_scheme_t  _scheme{CS_RK_NONE}; /*!< Selected RK scheme */
  char                    *_name{nullptr};      /*!< Associated equation's or
                                                     field name */
  int                      _n_stages{0};        /*!< Number of stages */
  int                      _i_stage{0};         /*!< Current stage index */
  cs_span<cs_real_t>       _dt;                 /*!< View of time step array */

  int                      _stride{0};          /*!< Variable stride */
  cs_lnum_t                _n_elts{0};          /*!< Number of computational
                                                     elements */

  /* variable storage */
  cs_array_2d<cs_real_t>   _u_old;              /*!< Variable at beginning of
                                                     the time step */
  cs_array_2d<cs_real_t>   _u_new;              /*!< Updated variable */
  cs_array_3d<cs_real_t>   _rhs_stages;         /*!< RHS temporary storage */
  cs_array<cs_real_t>      _mass;               /*!< Mass array for the
                                                     variables' resolution. */

  /* Coeffients adapted from the Butcher tableau */
  cs_array_2d<cs_real_t>   _a;            /*!< a RK coefficient */
  cs_array<cs_real_t>      _c;            /*!< c RK coefficient */
};
}

using cs_runge_kutta_integrator_t = cs::runge_kutta_integrator;
#endif

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
 * \brief Perform one Runge-Kutta staging for the potential in the NS system.
 *
 * \param[in]      ctx         Reference to dispatch context
 * \param[in,out]  rk          pointer to a Runge-Kutta integrator
 * \param[in,out]  phi_stage   high level variable array along stages
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_staging_potential(cs_dispatch_context          &ctx,
                                 cs_runge_kutta_integrator_t  *rk,
                                 cs_real_t                    *phi_stage);

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
  // Sanity check
  assert(rk != nullptr);
  assert(stride == rk->stride());

  // get the current stage index
  const int i_stg = rk->i_stage();

  cs_real_t *rhs = rk->get_rhs_stage_sub_array(i_stg);

  /* Allocate non reconstructed face value only if presence of limiter */

  /* We compute the total explicit balance. */

  int  inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
  int  imasac = 0;

  eqp->theta = 1;

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
    return rk->is_active();
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

  return rk->is_staging();
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
