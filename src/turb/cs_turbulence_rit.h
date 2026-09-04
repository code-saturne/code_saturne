#ifndef CS_TURBULENCE_RIT_H
#define CS_TURBULENCE_RIT_H

/*============================================================================
 * Turbulence transport equation.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_math.h"
#include <cmath>

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the divergence of turbulent flux to a scalar transport equation.
 *
 * \param[in]   field_id  transported field id
 * \param[in]   xcpp      Cp
 * \param[out]  vistet    diffusivity tensor
 * \param[out]  rhs       right hand side to update
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rit_div(const int        field_id,
                      const cs_real_t  xcpp[],
                      cs_real_t        vistet[][6],
                      cs_real_t        rhs[]);

/*---------------------------------------------------------------------------*/
/* Portable safe math helpers in device code                                 */
/*---------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE cs_real_t
_safe_exp(cs_real_t x)
{
  if (x > 700.)
    return exp(700.);
  if (x < -700.)
    return 0.;
  return exp(x);
}

static inline CS_F_HOST_DEVICE cs_real_t
_expm1_ratio(cs_real_t x)
{
  if (fabs(x) < 1.e-6)
    return 1. + 0.5*x + (1./6.)*x*x;
  return expm1(x)/x;
}

static inline CS_F_HOST_DEVICE cs_real_t
_log1p_ratio(cs_real_t x)
{
  if (fabs(x) < 1.e-6)
    return 1. - 0.5*x + (1./3.)*x*x;
  return log1p(x)/x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief 8-point Gauss-Legendre quadrature on [0, t].
 */
/*----------------------------------------------------------------------------*/

template <typename F>
static inline CS_F_HOST_DEVICE cs_real_t
_rij_source_gauss8(cs_real_t t, F&& f)
{
  static const cs_real_t gx[8] = {
    -0.960289856497536, -0.796666477413627,
    -0.525532409916329, -0.183434642495650,
     0.183434642495650,  0.525532409916329,
     0.796666477413627,  0.960289856497536
  };

  static const cs_real_t gw[8] = {
    0.101228536290376, 0.222381034453374,
    0.313706645877887, 0.362683783378362,
    0.362683783378362, 0.313706645877887,
    0.222381034453374, 0.101228536290376
  };

  if (t <= 0.)
    return 0.;

  const cs_real_t half_t = 0.5*t;
  cs_real_t sum = 0.;

  for (int i = 0; i < 8; i++) {
    const cs_real_t s = half_t*(1. + gx[i]);
    const cs_real_t value = f(s);
    sum += gw[i]*value;
  }

  return half_t*sum;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prescribed time scale and dissipative clock.
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE void
_source_time_stepping_tau_chi(cs_real_t   t,
                              cs_real_t   a_co,
                              cs_real_t   gamma_0,
                              cs_real_t   tau0,
                              cs_real_t  *tau_t,
                              cs_real_t  *chi_t)
{
  if (fabs(gamma_0) < 1.e-14) {
    *tau_t = tau0 + a_co*t;

    if (fabs(a_co) < 1.e-14)
      *chi_t = t/tau0;
    else
      *chi_t = log1p(a_co*t/tau0)/a_co;
  }
  else {
    const cs_real_t gamma_t = gamma_0*t;
    const cs_real_t e_minus = _safe_exp(-gamma_t);

    *tau_t = e_minus * (tau0 + a_co*t*_expm1_ratio(gamma_t));

    const cs_real_t z = t * _expm1_ratio(gamma_t);
    *chi_t = z * _log1p_ratio(a_co*z);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief J(dt; rate) = integral_0^dt exp(rate*chi(s)) ds.
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE cs_real_t
_source_time_stepping_ja(cs_real_t dt,
                         cs_real_t tau0,
                         cs_real_t a_co,
                         cs_real_t gamma_0,
                         cs_real_t rate)
{
  if (dt <= 0.)
    return 0.;

  if (fabs(a_co) < 1.e-12 && fabs(gamma_0) < 1.e-12) {
    if (fabs(rate) < 1.e-12)
      return dt;
    const cs_real_t gamma_rate = rate*dt/tau0;
    return dt*_expm1_ratio(gamma_rate);
  }

  return _rij_source_gauss8(
    dt,
    [=] CS_F_HOST_DEVICE (cs_real_t s) -> cs_real_t {
      cs_real_t tau_s;
      cs_real_t chi_s;
      _source_time_stepping_tau_chi(s, a_co, gamma_0, tau0, &tau_s, &chi_s);
      return _safe_exp(rate*chi_s);
    });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Derivative of J with respect to its rate.
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE cs_real_t
_source_time_stepping_dja_drate(cs_real_t dt,
                                cs_real_t tau0,
                                cs_real_t a_co,
                                cs_real_t gamma_0,
                                cs_real_t rate)
{
  if (dt <= 0.)
    return 0.;

  return _rij_source_gauss8(
    dt,
    [=] CS_F_HOST_DEVICE (cs_real_t s) -> cs_real_t {
      cs_real_t tau_s;
      cs_real_t chi_s;
      _source_time_stepping_tau_chi(s, a_co, gamma_0, tau0, &tau_s, &chi_s);
      const cs_real_t exponent = rate*chi_s;
      return chi_s*_safe_exp(exponent);
    });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief F_b(dt) = J(dt; (C_R - C_theta)/2).
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE cs_real_t
_source_time_stepping_Fb(cs_real_t dt,
                         cs_real_t tau0,
                         cs_real_t a_co,
                         cs_real_t gamma_0,
                         cs_real_t cr,
                         cs_real_t ctheta)
{
  return _source_time_stepping_ja(dt,
                                  tau0,
                                  a_co,
                                  gamma_0,
                                  0.5*(cr - ctheta));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief L_k(dt) = J(dt; 1 - (C_R + C_theta)/2).
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE cs_real_t
_source_time_stepping_Lk(cs_real_t dt,
                         cs_real_t tau0,
                         cs_real_t a_co,
                         cs_real_t gamma_0,
                         cs_real_t cr,
                         cs_real_t ctheta)
{
  return _source_time_stepping_ja(dt,
                                  tau0,
                                  a_co,
                                  gamma_0,
                                  1. - 0.5*(cr + ctheta));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief H_k(dt).
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE cs_real_t
_source_time_stepping_hk(cs_real_t dt,
                         cs_real_t tau0,
                         cs_real_t a_co,
                         cs_real_t gamma_0,
                         cs_real_t cr,
                         cs_real_t ctheta)
{
  const cs_real_t buoyancy_diff = 0.5*(cr - ctheta);
  const cs_real_t kinetic_rate  = 1. - 0.5*(cr + ctheta);
  const cs_real_t gamma_scale   = fmax(1., fabs(gamma_0));

  if (fabs(gamma_0) <= 1.e-12*gamma_scale) {
    const cs_real_t denom = a_co + buoyancy_diff;
    const cs_real_t scale =
      fmax(1., fmax(fabs(a_co), fabs(buoyancy_diff)));

    if (fabs(denom) > 1.e-10*scale) {
      const cs_real_t upper =
        _source_time_stepping_ja(dt,
                                 tau0,
                                 a_co,
                                 gamma_0,
                                 kinetic_rate + denom);

      const cs_real_t lower =
        _source_time_stepping_ja(dt,
                                 tau0,
                                 a_co,
                                 gamma_0,
                                 kinetic_rate);

      return tau0/denom*(upper - lower);
    }

    return _source_time_stepping_dja_drate(dt,
                                           tau0,
                                           a_co,
                                           gamma_0,
                                           kinetic_rate);
  }

  return _rij_source_gauss8(
    dt,
    [=] CS_F_HOST_DEVICE (cs_real_t s) -> cs_real_t {
      cs_real_t tau_s;
      cs_real_t chi_s;

      _source_time_stepping_tau_chi(s,
                                    a_co,
                                    gamma_0,
                                    tau0,
                                    &tau_s,
                                    &chi_s);

      const cs_real_t fb_s =
        _source_time_stepping_Fb(s,
                                 tau0,
                                 a_co,
                                 gamma_0,
                                 cr,
                                 ctheta);

      return _safe_exp(kinetic_rate*chi_s)*fb_s;
    });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Common source-step engine for frozen and prescribed variable tau.
 */
/*----------------------------------------------------------------------------*/

static inline CS_F_HOST_DEVICE void
_source_time_stepping_engine(cs_real_t          cr,
                             cs_real_t          ctheta,
                             cs_real_t          ceps2,
                             cs_real_t          a_co,
                             cs_real_t          gamma_0,
                             bool               variable_tau_mode,
                             cs_real_t          beta,
                             const cs_real_t    grav[3],
                             cs_real_t          dt,
                             const cs_real_6_t  r0,
                             const cs_real_3_t  qtheta0,
                             cs_real_t          theta2_0,
                             cs_real_t          eps0,
                             cs_real_6_t        r1,
                             cs_real_3_t        qtheta1,
                             cs_real_t         *theta2_1,
                             cs_real_t         *eps1)
{
  const cs_real_t tau0 = 1.0/eps0;

  cs_real_t tau_dt;
  cs_real_t chi_dt;

  if (!variable_tau_mode) {
    tau_dt = tau0;
    chi_dt = dt/tau0;
  }
  else {
    _source_time_stepping_tau_chi(dt,
                                  a_co,
                                  gamma_0,
                                  tau0,
                                  &tau_dt,
                                  &chi_dt);
  }

  *eps1 = 1.0/tau_dt;

  const cs_real_t fb_dt =
    _source_time_stepping_Fb(dt,
                             tau0,
                             a_co,
                             gamma_0,
                             cr,
                             ctheta);

  const cs_real_t lk_dt =
    _source_time_stepping_Lk(dt,
                             tau0,
                             a_co,
                             gamma_0,
                             cr,
                             ctheta);

  const cs_real_t hk_dt =
    _source_time_stepping_hk(dt,
                             tau0,
                             a_co,
                             gamma_0,
                             cr,
                             ctheta);

  const cs_real_t buoyancy = cs_math_3_dot_product(qtheta0, grav);
  const cs_real_t b_exp = _safe_exp((1. - cr)*chi_dt);

  *theta2_1 = b_exp * (theta2_0 + 2.0*buoyancy*fb_dt);

  const cs_real_t g_term = lk_dt*theta2_0 + 2.0*buoyancy*hk_dt;

  for (cs_lnum_t ij = 0; ij < 3; ij++) {
    qtheta1[ij] = _safe_exp(-0.5*(cr + ctheta)*chi_dt)
                  * (qtheta0[ij] - grav[ij]*beta*g_term);
  }

  const cs_real_t exp_cr = _safe_exp(-cr*chi_dt);
  const cs_real_t exp_eps = _safe_exp(-ceps2*chi_dt);

  cs_real_t g_delta = 0.;
  if (variable_tau_mode) {
    const cs_real_t ceps3 = 1. - gamma_0*tau0/((1. - ceps2)*beta*buoyancy);
    g_delta = (1. - ceps3)*beta;
  }

  const cs_real_t j_ceps2 =
    _source_time_stepping_ja(dt, tau0, a_co, gamma_0, 1. - ceps2);

  const cs_real_t j_cr =
    _source_time_stepping_ja(dt, tau0, a_co, gamma_0, 1. - cr);

  const cs_real_t h_ceps2 =
    _source_time_stepping_hk(dt, tau0, a_co, gamma_0, cr, 2.0*ceps2 - cr);

  const cs_real_t h_cr =
    _source_time_stepping_hk(dt, tau0, a_co, gamma_0, cr, cr);

  for (cs_lnum_t ij = 0; ij < 6; ij++) {
    cs_real_t source_term = 0.;
    if (ij < 3) {
      source_term =
        - (2.0/3.0) * (1.0 - exp_eps) * (ij == 0 ? r0[0]+r0[1]+r0[2] : 0.)
        - 2.0*grav[ij]*beta*qtheta0[ij]*j_cr
        + 2.0*cs_math_pow2(grav[ij])*beta*g_term*h_cr;
    }
    else {
      cs_lnum_t i = (ij == 3) ? 0 : (ij == 4 ? 1 : 2);
      cs_lnum_t j = (ij == 3) ? 1 : (ij == 4 ? 2 : 0);
      source_term =
        - beta*(grav[i]*qtheta0[j] + grav[j]*qtheta0[i])*j_cr
        + 2.0*grav[i]*grav[j]*beta*g_term*h_cr;
    }

    cs_real_t source_eps = 0.;
    if (ij < 3) {
      source_eps =
        - (2.0/3.0) * g_delta * buoyancy * j_ceps2 * (ij == 0 ? 1. : 0.)
        - 2.0 * g_delta * grav[ij] * qtheta0[ij] * h_ceps2;
    }
    else {
      cs_lnum_t i = (ij == 3) ? 0 : (ij == 4 ? 1 : 2);
      cs_lnum_t j = (ij == 3) ? 1 : (ij == 4 ? 2 : 0);
      source_eps =
        - g_delta * (grav[i]*qtheta0[j] + grav[j]*qtheta0[i]) * h_ceps2;
    }

    r1[ij] = exp_cr*r0[ij] + source_term + source_eps;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief CS_TURB_RIJ_SOURCE_TS_EXPONENTIAL: exact frozen-tau integration.
 */
/*----------------------------------------------------------------------------*/

inline CS_F_HOST_DEVICE void
cs_turbulence_rit_source_step_frozen_tau(cs_real_t          cr,
                                         cs_real_t          ctheta,
                                         cs_real_t          ceps2,
                                         cs_real_t          beta,
                                         const cs_real_t    grav[3],
                                         cs_real_t          dt,
                                         const cs_real_6_t  r0,
                                         const cs_real_3_t  qtheta0,
                                         cs_real_t          theta2_0,
                                         cs_real_t          eps0,
                                         cs_real_6_t        r1,
                                         cs_real_3_t        qtheta1,
                                         cs_real_t         *theta2_1,
                                         cs_real_t         *eps1)
{
  _source_time_stepping_engine(cr,
                               ctheta,
                               ceps2,
                               0.,
                               0.,
                               false,
                               beta,
                               grav,
                               dt,
                               r0,
                               qtheta0,
                               theta2_0,
                               eps0,
                               r1,
                               qtheta1,
                               theta2_1,
                               eps1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief CS_TURB_RIJ_SOURCE_TS_VAR_TAU: variable-tau trajectory integration.
 */
/*----------------------------------------------------------------------------*/

inline CS_F_HOST_DEVICE void
cs_turbulence_rit_source_step_variable_tau(cs_real_t          cr,
                                           cs_real_t          ctheta,
                                           cs_real_t          ceps2,
                                           cs_real_t          ceps3,
                                           cs_real_t          beta,
                                           const cs_real_t    grav[3],
                                           cs_real_t          dt,
                                           const cs_real_6_t  r0,
                                           const cs_real_3_t  qtheta0,
                                           cs_real_t          theta2_0,
                                           cs_real_t          eps0,
                                           cs_real_6_t        r1,
                                           cs_real_3_t        qtheta1,
                                           cs_real_t         *theta2_1,
                                           cs_real_t         *eps1)
{
  const cs_real_t k0 =
    0.5*(r0[0] + r0[1] + r0[2]);

  const cs_real_t qg0 =
    cs_math_3_dot_product(qtheta0, grav);

  const cs_real_t a_co = ceps2 - 1.;

  const cs_real_t gamma_0 =
    (1. - ceps3)*beta*qg0/k0;

  _source_time_stepping_engine(cr,
                               ctheta,
                               ceps2,
                               a_co,
                               gamma_0,
                               true,
                               beta,
                               grav,
                               dt,
                               r0,
                               qtheta0,
                               theta2_0,
                               eps0,
                               r1,
                               qtheta1,
                               theta2_1,
                               eps1);
}

#endif /* CS_TURBULENCE_RIT_H */
