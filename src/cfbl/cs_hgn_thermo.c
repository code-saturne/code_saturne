/*============================================================================
 * Thermodynamic for the compressible homogeneous two-phase flow model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_math.h"
#include "cs_hgn_phase_thermo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hgn_thermo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_hgn_thermo.c
 *
 *  \brief Thermodynamic of a compressible homogeneous two-phase flow
 *
 *  In order to define the thermodynamical behaviour of the compressible
 *  homogeneous two-phase flow, an Equation Of States (EOS) has to be specified
 *  for both phases (see the file \ref cs_hgn_phase_thermo.c).
 *
 *  In the present file can be found the computation of the mixture properties:
 *    - pressure
 *    - temperature
 *    - sound speed in the mixture
 *    - specific entropy
 *    - specific internal energy.
 *
 *  These properties are built on the basis of the assumptions:
 *    - that the extensive mixture entropy is equal to the sum of the two
 *      extensive phasic entropies
 *    - that the Gibbs relation holds within each single phase.
 *
 *  These assumptions then permit to compute a Gibbs relation for the
 *  mixture and thus they allow to define the pressure and temperature law
 *  for the mixture.
 */
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static local variables
 *============================================================================*/

/* limit values of fraction used to define single-phase and two-phase regimes */
static cs_real_t _eps_fraction_lim = 1.e-12;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * functional for computation of saturation temperature
 * f(T) = mu_1/T-mu-2/T
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_tsat_function(cs_real_t tp, cs_real_t pr, cs_real_t *pdfunc)
{
  cs_real_t tplus = 1.00001*tp;

  cs_real_t mu1 =  cs_hgn_phase_thermo_internal_energy_tp(tp, pr, 0)
                 + cs_hgn_phase_thermo_specific_volume_tp(tp, pr, 0)*pr
                 - tp * cs_hgn_phase_thermo_entropy_tp(tp, pr, 0);
  cs_real_t mu2 =  cs_hgn_phase_thermo_internal_energy_tp(tp, pr, 1)
                 + cs_hgn_phase_thermo_specific_volume_tp(tp, pr, 1)*pr
                 - tp * cs_hgn_phase_thermo_entropy_tp(tp, pr, 1);

  cs_real_t res = mu1 - mu2;

  mu1 =  cs_hgn_phase_thermo_internal_energy_tp(tplus, pr, 0)
       + cs_hgn_phase_thermo_specific_volume_tp(tplus, pr, 0)*pr
       - tplus * cs_hgn_phase_thermo_entropy_tp(tplus, pr, 0);
  mu2 =  cs_hgn_phase_thermo_internal_energy_tp(tplus, pr, 1)
       + cs_hgn_phase_thermo_specific_volume_tp(tplus, pr, 1)*pr
       - tplus * cs_hgn_phase_thermo_entropy_tp(tplus, pr, 1);

  cs_real_t dfunc = ((mu1-mu2)-res) / (tplus-tp);

  *pdfunc = dfunc;
  return res;
}

/*----------------------------------------------------------------------------
 * Function used to compute saturation temperature with a secant method
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_secant_tsat(cs_real_t pr, cs_real_t tini)
{
  cs_real_t dfn;
  cs_real_t tn = tini;
  cs_real_t fn = _tsat_function(tn, pr, &dfn);

  cs_real_t tnp = tini*1.0001;
  cs_real_t fnp = _tsat_function(tnp, pr, &dfn);
  dfn = (fnp - fn)/(tnp - tn);

  int iterm = 100;
  cs_real_t eps = 1.0e-10;
  for (int iter = 0; iter < iterm; iter++) {
    if (CS_ABS(fn) < eps) break;
    tnp = tnp - fnp / dfn;
    fnp = _tsat_function(tnp, pr, &dfn);
    dfn = (fnp - fn) / (tnp - tn);
    tn=tnp;
    fn=fnp;
  }

  return tn;
}

/*----------------------------------------------------------------------------
 * Computation of the energy fraction, the mixture pressure and specific
 * energy  with respect to the specific entropy, the specific volume, the mass,
 * volume and entropy fractions.
 *
 * Used in the computation of the sound speed in the mixture.
 *
 * alpha   --> volume fraction
 * y       --> mass fraction
 * beta    --> entropy fraction
 * s       --> specific entropy
 * v       --> specific volume
 * pz     <--  energy fraction
 * pe     <--  specific energy
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_mix_pressure_sv(cs_real_t alpha,
                 cs_real_t y,
                 cs_real_t beta,
                 cs_real_t s,
                 cs_real_t v,
                 cs_real_t *pz,
                 cs_real_t *pe)
{
  cs_real_t eps = _eps_fraction_lim;
  cs_real_t pr, z, e;

  if (v < 0.)
    bft_error(__FILE__, __LINE__, 0,
              _("Input of mix pressure computation with respect to specific "
                "entropy and specific volume:\n mix specific volume v < 0\n"));

  /* single-phase : phase 2 */
  if (alpha < eps || y < eps || beta < eps) {

    cs_real_t s_2 = s;
    cs_real_t v_2 = v;

    cs_real_t e_2 = cs_hgn_phase_thermo_internal_energy_sv(s_2, v_2, 1);

    z = y;
    e = e_2;
    pr = cs_hgn_phase_thermo_pressure_ve(v_2, e_2, 1);

  /* single-phase : phase 1 */
  } else if (   1. - alpha < eps || 1 - y < eps || 1. - beta < eps) {

    cs_real_t s_1 = s;
    cs_real_t v_1 = v;

    cs_real_t e_1 = cs_hgn_phase_thermo_internal_energy_sv(s_1, v_1, 0);

    z = y;
    e = e_1;
    pr = cs_hgn_phase_thermo_pressure_ve(v_1, e_1, 0);

  /* two-phase mixture */
  } else {

    cs_real_t v_1 = alpha * v/y;
    cs_real_t v_2 = (1. - alpha) * v / (1. - y);
    cs_real_t s_1 = beta * s / y;
    cs_real_t s_2 = (1. - beta) * s / (1. - y);

    cs_real_t e_1 = cs_hgn_phase_thermo_internal_energy_sv(s_1, v_1, 0);
    cs_real_t e_2 = cs_hgn_phase_thermo_internal_energy_sv(s_2, v_2, 1);

    e = y * e_1 + (1. - y) * e_2;

    if (e < 0.)
      bft_error(__FILE__, __LINE__, 0,
                _("While computing mix pressure with respect to specific "
                  "entropy and specific volume:\n mix internal energy e < 0\n"));

    z = y * e_1 / e;
    cs_real_t tp1 = cs_hgn_phase_thermo_temperature_ve(v_1, e_1, 0);
    cs_real_t tp2 = cs_hgn_phase_thermo_temperature_ve(v_2, e_2, 1);
    cs_real_t invt = z / tp1 + (1. - z) / tp2;

    if (isnan(invt)) {
      bft_printf(_("In _mix_pressure_sv : 1/temperature NAN\n"));
    }

    cs_real_t tp = 1. / invt;

    if (tp < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("While computing mix pressure with respect to specific "
                  "entropy and specific volume:\n mix temperature T < 0\n"));

    cs_real_t pdivt =  alpha      * cs_hgn_phase_thermo_pressure_ve(v_1, e_1, 0) / tp1
                     + (1.-alpha) * cs_hgn_phase_thermo_pressure_ve(v_2, e_2, 1) / tp2;
    pr = pdivt * tp;
  }

  *pe = e;
  *pz = z;

  return pr;
}

/*----------------------------------------------------------------------------
 * Function \f$F\f$ used for the computation of the equilibrium fractions
 *
 * Two different (and equivalent) forms of the functions are available
 * through the use of the parameter "form".
 * The two forms of the functions aim at dealing with possible
 * ill-conditioning issues:
 * - form 1 for two-phase non evanescent
 * - form 2 for two-phase evanescent
 *
 * \param[in]  e      specific internal energy
 * \param[in]  v      specific volume
 * \param[in]  p      pressure
 * \param[in]  form   form parameter
 *
 * \return \f$ F(e,v,p) \f$
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_eq_function(cs_real_t e,
             cs_real_t v,
             cs_real_t p,
             int       form)
{
  cs_real_t res1, res2;
  cs_real_t tsat = cs_hgn_thermo_saturation_temp(p);

  if (form == 1) {

    res1 =  (e - cs_hgn_phase_thermo_internal_energy_tp(tsat, p, 1))
          * (  cs_hgn_phase_thermo_specific_volume_tp(tsat, p, 0)
             - cs_hgn_phase_thermo_specific_volume_tp(tsat, p, 1));
    res2 =  (v - cs_hgn_phase_thermo_specific_volume_tp(tsat, p, 1))
          * (  cs_hgn_phase_thermo_internal_energy_tp(tsat, p, 0)
             - cs_hgn_phase_thermo_internal_energy_tp(tsat, p, 1));

  } else if (form == 2) {

    res1 =  (e - cs_hgn_phase_thermo_internal_energy_tp(tsat, p, 1))
          / (  cs_hgn_phase_thermo_internal_energy_tp(tsat, p, 0)
             - cs_hgn_phase_thermo_internal_energy_tp(tsat, p, 1));
    res2 =  (v - cs_hgn_phase_thermo_specific_volume_tp(tsat, p, 1))
          / (  cs_hgn_phase_thermo_specific_volume_tp(tsat, p, 0)
             - cs_hgn_phase_thermo_specific_volume_tp(tsat, p, 1));

  } else {

    bft_error(__FILE__, __LINE__, 0,
              _("Unknown form for equilibrium function.\n"));

  }

  return res1-res2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dichotomy algorithm for the computation of the equilibrium fractions.
 *
 * It is based on a Dichotomy algorithm on the equilibrium function (which can
 * have several forms) to search (if it exists) a solution between the
 * pressures \f$p_a\f$ and \f$p_b\f$ (\f$ pa < pb \f$).
 *
 * \param[in]  e           specific internal energy
 * \param[in]  v           specific volume
 * \param[in]  pa          minimum pressure (Pa)
 * \param[in]  pb          maximum pressure (Pa)
 * \param[out] palpha_eq   equilibrium volume fraction
 * \param[out] py_eq       equilibrium mass fraction
 * \param[out] pz_eq       equilibrium energy fraction
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_dicho_eq(cs_real_t  e,
          cs_real_t  v,
          cs_real_t  pa,
          cs_real_t  pb,
          cs_real_t *palpha_eq,
          cs_real_t *py_eq,
          cs_real_t *pz_eq)
{
  cs_real_t alpha_eq, y_eq, z_eq;

  cs_real_t pmin = pa;
  cs_real_t pmax = pb;

  int form = 1;
  cs_real_t fmin = _eq_function(e, v, pmin, form);
  cs_real_t fmax = _eq_function(e, v, pmax, form);

  cs_real_t ptmp;

  /* no solution */
  if (fmin * fmax > 0)
    form = -1;

  /* start two-phase dichotomy search */
  if (form > 0) {

    cs_real_t ftmp, ptmpm1;
    int imax = 100;

    for (int iter = 0; iter <= imax; iter++) {

      if (iter >= 1) ptmpm1 = ftmp;

      ptmp = 0.5*(pmin+pmax);
      ftmp = _eq_function(e, v, ptmp, form);

      /* FIXME avoid hard coded precision */
      if (iter >= 1 && CS_ABS(ptmp-ptmpm1) < 1.e-8*CS_ABS(ptmp)) break;
      if (CS_ABS(ftmp) < 1.e-8) break;

      if (fmin * ftmp < 0) {
        pmax = ptmp;
        fmax = ftmp;
      }
      else if (fmax * ftmp <= 0) {
        pmin = ptmp;
        fmin = ftmp;
      }
      else {
        bft_error(__FILE__, __LINE__, 0,
                  _("While performing dichotomy search on equilibrium"
                    " function\n"));
      }
    }

    cs_real_t ts = cs_hgn_thermo_saturation_temp(ptmp);

    y_eq =  (v - cs_hgn_phase_thermo_specific_volume_tp(ts, ptmp, 1))
          / (  cs_hgn_phase_thermo_specific_volume_tp(ts, ptmp, 0)
             - cs_hgn_phase_thermo_specific_volume_tp(ts, ptmp, 1));

    alpha_eq = y_eq*cs_hgn_phase_thermo_specific_volume_tp(ts, ptmp, 0) / v;
    z_eq = y_eq*cs_hgn_phase_thermo_internal_energy_tp(ts, ptmp, 0) / e;

  } else {

    alpha_eq = -1.;
    y_eq = -1.;
    z_eq = -1.;
    ptmp = -1.;

  }

  *palpha_eq = alpha_eq;
  *py_eq = y_eq;
  *pz_eq = z_eq;

  return ptmp;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the temperature at saturation with respect to the
 *        pressure.
 *
 * Compute the temperature at saturation \f$T_{sat}\f$ with respect to the
 * pressure \f$P\f$. It corresponds to the temperature for which
 * the pressure, temperature and chemical potential are equal in both phases.
 * It reduces to solve the equality of the phasic chemical potential:
 * \f[
 *   P \rightarrow T_{sat} \/ \mu_1(T_{sat},P)=\mu_2(T_{sat},P)
 * \f]
 * This equality is solved using a Newton method.
 *
 * \param[in]     pr       pressure
 *
 * \return temperature at saturation
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_thermo_saturation_temp(cs_real_t pr)
{
  /* Initial guesses for the Newton algorithm.
   * tab_tsat contains 101 pressure values from 700 Pa to 25000000 Pa with an
   * increment of 249993. This small lookup table (piecewise constant) is then
   * used as an initial guess for the Newton algorithm that compute Tsat.*/
  cs_real_t tab_tsat[101] = {
    275.040788823, 400.672644776, 425.0456133,440.945408602, 453.063793165,
    462.992518065, 471.461895813, 478.896713175, 485.541070517, 491.572761469,
    497.106364218, 502.228041699, 507.002478334, 511.480126122, 515.700428738,
    519.695831821, 523.492159312, 527.109853231, 530.568952531, 533.884656876,
    537.068860536, 540.13333692, 543.089241621, 545.944512744, 548.706058795,
    551.380794855, 553.975546325, 556.495978335, 558.946076123, 561.329785387,
    563.651051827, 565.913821138, 568.122000266, 570.278579733, 572.385768942,
    574.445758103, 576.460737428, 578.432897126, 580.364427409, 582.257571137,
    584.114077596, 585.935272099, 587.72238163, 589.476633171, 591.199253705,
    592.891470214, 594.554509681, 596.189599087, 597.797919408, 599.38022613,
    600.937267041, 602.470053718, 603.979298284, 605.465639951, 606.929717933,
    608.372171441, 609.79363969, 611.194761892, 612.576177259, 613.938544532,
    615.282247716, 616.607755338, 617.915532784, 619.206042845, 620.479714099,
    621.736925444, 622.978077771, 624.20354135, 625.413681693, 626.608824263,
    627.789349263, 628.955496708, 630.10762033, 631.245999792, 632.370928664,
    633.482659645, 634.581450284, 635.667558128, 636.741240647, 637.802733775,
    638.852255407, 639.890023444, 640.916255788, 641.931160006, 642.934941835,
    643.927787539, 644.909898305, 645.881429323, 646.842724477, 647.783039797,
    648.710037116, 649.63489674, 650.557302478, 651.476996344, 652.393682757,
    653.307179303, 654.217316457, 655.123897609, 656.026880438, 656.926042113,
    657.821270598 };

  cs_real_t dp = 249993.;
  cs_real_t pmin = 700.;
  int ip = (int)((pr-pmin)/dp);

  if (ip < 0) ip = 0;
  if (ip > 100) ip = 100;

  /* Call to the Newton algorithm with tabTsat[ip] as an initial guess.*/
  /* FIXME It should be replaced by a value computed from tabTsat[] as a piecewise
     linear approximation. It's easy and will provide a more precise initial guess
     and it will thus lower the number of iterations in the Newton algorithm. */
  cs_real_t tsat = _secant_tsat(pr, tab_tsat[ip]);

  if (isnan(tsat))
    tsat = -1.;

  return tsat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation mixture pressure and temperature from volume, mass,
 *        energy fractions, as well as specific energy and specific volume.
 *
 * Following relations are used, that rely on phasic pressures and temperatures:
 * \f[
 *   \dfrac{1}{T} = \dfrac{\dd s}{\dd e}_{|\tau,\alpha,y,z}
 * \f]
 * \f[
 *   \dfrac{P}{T} = \dfrac{\dd s}{\dd \tau}_{|e,\alpha,y,z}
 * \f]
 *
 * \param[in]     alpha   volume fraction
 * \param[in]     y       mass fraction
 * \param[in]     z       energy fraction
 * \param[in]     e       specific energy
 * \param[in]     v       specific volume
 * \param[out]    ptp     pointer to mixture temperature
 * \param[out]    ppr     pointer to mixture pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_hgn_thermo_pt(cs_real_t  alpha,
                 cs_real_t  y,
                 cs_real_t  z,
                 cs_real_t  e,
                 cs_real_t  v,
                 cs_real_t *ptp,
                 cs_real_t *ppr)
{
  cs_real_t tp, pr;

  if (v <= 0.)
    bft_error(__FILE__, __LINE__, 0,
              _("Input of mix pressure and temperature computation with "
                "respect to specific energy and specific volume:\n"
                "specific volume <= 0\n"));

  if (e <= 0.)
    bft_error(__FILE__, __LINE__, 0,
              _("Input of mix pressure and temperature computation with "
                "respect to specific energy and specific volume:\n"
                "specific energy <= 0\n"));

  /* single-phase : phase 2 */
  if (y < _eps_fraction_lim || z < _eps_fraction_lim) {

    tp = cs_hgn_phase_thermo_temperature_ve(v, e, 1);
    if (tp < 0.)
      bft_error(__FILE__, __LINE__, 0,
                _("Single-phase regime - phase 2: temperature < 0\n"));
    pr = cs_hgn_phase_thermo_pressure_ve(v, e, 1);

  /* single-phase : phase 1 */
  } else if (1-y < _eps_fraction_lim || 1-z < _eps_fraction_lim) {

    tp = cs_hgn_phase_thermo_temperature_ve(v, e, 0);
    if (tp < 0.)
      bft_error(__FILE__, __LINE__, 0,
                _("Single-phase regime - phase 1: temperature < 0\n"));
    pr = cs_hgn_phase_thermo_pressure_ve(v, e, 0);

  } else {

    cs_real_t e1 = z*e/y;
    cs_real_t v1 = alpha*v/y;
    cs_real_t e2 = (1.-z)*e/(1.-y);
    cs_real_t v2 = (1.-alpha)*v/(1.-y);

    cs_real_t tp1 = cs_hgn_phase_thermo_temperature_ve(v1, e1, 0);
    cs_real_t tp2 = cs_hgn_phase_thermo_temperature_ve(v2, e2, 1);
    cs_real_t p1 = cs_hgn_phase_thermo_pressure_ve(v1, e1, 0);
    cs_real_t p2 = cs_hgn_phase_thermo_pressure_ve(v2, e2, 1);

    cs_real_t invt = z/tp1 + (1.-z)/tp2;
    if (isnan(invt))
      bft_printf(_("cs_hgn_thermo_pt() : 1.0/temperature NAN  (two-phase)\n"));

    tp = 1./invt;

    if (tp < 0.)
      bft_error(__FILE__, __LINE__, 0,
                _("Two-phase regime: mixture temperature < 0\n"));

    cs_real_t pdivt = alpha*p1/tp1+(1.-alpha)*p2/tp2;
    pr = pdivt*tp;

  }

  if (isnan(tp))
    bft_printf(_("cs_hgn_thermo_pt() : temperature NAN\n"));

  if (isnan(pr))
    bft_printf(_("cs_hgn_thermo_pt() : pressure NAN\n"));

  *ppr = pr;
  *ptp = tp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the square of the sound speed in the mixture.
 *
 * The sound speed may be computed through the Hessian matrices of the specific
 * phasic entropies in the plane (\f$\tau\f$, \f$e\f$).
 * \f$\tau\f$ stands for specific volume, and \f$e\f$ for specific energy.
 * The sound speed is here estimated using the plane (\f$\tau\f$, \f$s\f$).
 * \f$s\f$ stands for specific entropy.
 * We use the definition
 * \f\[
 *   c^2 = -\tau^2 \der{P}{\tau}_{|s,y}
 * \f\].
 * This relation is computed by a finite difference.
 *
 * \param[in]     alpha     volume fraction
 * \param[in]     y         mass fraction
 * \param[in]     z         energy fraction
 * \param[in]     P         pressure
 * \param[in]     v         specific volume
 *
 * \return square of the sound speed.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_thermo_c2(cs_real_t alpha,
                 cs_real_t y,
                 cs_real_t z,
                 cs_real_t P,
                 cs_real_t v)
{
  cs_real_t e = cs_hgn_thermo_ie(alpha, y, z, P, v);

  cs_real_t s, beta;

  /* single-phase : phase 2 */
  if (y <= _eps_fraction_lim) {

    s = cs_hgn_phase_thermo_entropy_ve(v, e, 1);
    beta = y;

  }
  /* single-phase : phase 1 */
  else if (1.-y <= _eps_fraction_lim) {

    s = cs_hgn_phase_thermo_entropy_ve(v, e, 0);
    beta = y;

  }
  /* two-phase */
  else {

    cs_real_t s1 = cs_hgn_phase_thermo_entropy_ve(alpha*v/y, z*e/y, 0);
    cs_real_t s2 = cs_hgn_phase_thermo_entropy_ve((1.-alpha)*v/(1.-y), (1.-z)*e/(1.-y), 1);
    s = (y*s1+(1.-y)*s2);
    beta = y*s1/s;

  }

  /* square of sound celerity */
  cs_real_t dv = 1.e-3*v;
  cs_real_t celer = -v*v*( _mix_pressure_sv(alpha, y, beta, s, v+dv, &z, &e)
                          -_mix_pressure_sv(alpha, y, beta, s, v   , &z, &e))
                     / dv;

  if (isnan(celer))
    bft_printf(_("cs_hgn_thermo_c2() : NAN\n"));

  /* hyperbolicity test */
  if (celer < cs_math_epzero)
    bft_error(__FILE__, __LINE__, 0,
              _("Negative sound speed - hyperbolicity problem\n"));

  return celer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the specific internal energy with respect to the
 *        volume (\f$\alpha\f$), mass (\f$y\f$) and energy (\f$\z\f$) fractions,
 *        as well as the pressure and the specific volume \f$\tau\f$.
 *
 * It uses a quasi-Newton method to solve:
 * \f[
 *   \mcal{P}(\alpha, y, z, e, \tau) = P
 * \f]
 *
 * \param[in]     alpha     the volume fraction
 * \param[in]     y         the mass fraction
 * \param[in]     z         the energy fraction
 * \param[in]     pr         the pressure
 * \param[in]     v         the specific volume
 *
 * \return specific internal energy.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_thermo_ie(cs_real_t alpha,
                 cs_real_t y,
                 cs_real_t z,
                 cs_real_t pr,
                 cs_real_t v)
{
  /* pctau is the critical pressure, first intersection point of
     the curves tau_sat_k(P). It also corresponds with the point where Tsat(P)
     starts decreasing -> FIXME WARNING: it has to be changed with the thermo */
  cs_real_t pctau = 1.5665e8;

  /* Initial guess for the quasi-Newton method. */

  cs_real_t tsat = cs_hgn_thermo_saturation_temp(pctau);
  cs_real_t en = CS_MAX(cs_hgn_phase_thermo_internal_energy_tp(tsat, pctau, 0),
                        cs_hgn_phase_thermo_internal_energy_tp(tsat, pctau, 1));

  cs_real_t tp, pn;
  cs_hgn_thermo_pt(alpha, y, z, en, v, &tp, &pn);

  cs_real_t de = 1.e-2*en;

  /* quasi-Newton method iterations. */
  int imax = 1000;
  for (int iter = 0; iter < imax; ++iter) {
    cs_real_t fn = (pn-pr);

    /* FIXME avoid hard coded precision */
    if (CS_ABS(fn / pr) < 1.e-10)
      break;

    cs_real_t pnpde;
    cs_hgn_thermo_pt(alpha, y, z, en+de, v, &tp, &pnpde);
    cs_real_t df = (pnpde - pn) / de;

    /* FIXME avoid hard coded precision */
    if (CS_ABS(df) < 1.e-8)
      break;

    en = en-fn/df;
    cs_hgn_thermo_pt(alpha, y, z, en, v, &tp, &pn);
  }

  if (en < 0.)
    bft_error(__FILE__, __LINE__, 0,
              _("Negative specific internal energy e < 0\n"));

  return en;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the equilibrium fractions.
 *
 * The equilibrium fractions correspond to the definition of the mixture for
 * which one gets the pressure, temperature and chemical potential equilibrium.
 *
 * They are computed by using a Dichotomy algorithm on the function
 * characterizing the equilibrium (two forms available).
 *
 * The search for the equilibrium point is done in plane (P,T). Dichotomy is
 * performed on the pressure along the saturation curve.
 *
 * \param[in]  e          specific internal energy
 * \param[in]  v          specific volume
 * \param[out] palpha_eq  pointer to equilibrium volume fraction
 * \param[out] py_eq      pointer to equilibrium mass fraction
 * \param[out] pz_eq      pointer to equilibrium energy fraction
 */
/*----------------------------------------------------------------------------*/

void
cs_hgn_thermo_eq(cs_real_t  e,
                 cs_real_t  v,
                 cs_real_t *palpha_eq,
                 cs_real_t *py_eq,
                 cs_real_t *pz_eq)
{
  /* pctau is the critical pressure, first intersection point of
     the curves tau_sat_k(P). It also corresponds with the point where Tsat(P)
     starts decreasing -> FIXME WARNING: it has to be changed with the thermo */
  cs_real_t pctau = 1.5665e8;

  cs_real_t alpha_eq, y_eq, z_eq;

  cs_real_t pmin = 1.;

  /* search root in [pmin,0.5(pmin+pctau)] */
  cs_real_t pa = pmin;
  cs_real_t pb = 0.5*(pmin+pctau);

  cs_real_t alpha1, y1, z1;
  cs_real_t p1 = _dicho_eq(e, v, pa, pb, &alpha1, &y1, &z1);

  /* root not found in [pmin,0.5(pmin+pctau)] */
  if (   (alpha1 < 0. || alpha1 > 1.)
      || (y1     < 0. || y1     > 1.)
      || (z1     < 0. || z1     > 1.)) {

    /* search root in [0.5(pmin+pctau),pctau] */
    pa = 0.5*(pmin+pctau);
    pb = pctau;
    cs_real_t alpha2, y2, z2;
    cs_real_t p2 = _dicho_eq(e, v, pa, pb, &alpha2, &y2, &z2);

    /* root not found neither in [0.5(pmin+pctau),pctau] */
    if (   (alpha2 < 0. || alpha2 > 1.)
        || (y2     < 0. || y2     > 1.)
        || (z2     < 0. || z2     > 1.)) {

      /* single phase regime */

      /* pick a single-phase configuration */

      cs_real_t s_1 = cs_hgn_phase_thermo_entropy_ve(v, e, 0);
      cs_real_t s_2 = cs_hgn_phase_thermo_entropy_ve(v, e, 1);

      if (s_1 > s_2) {
        alpha_eq = 1.;
        y_eq = 1.;
        z_eq = 1.;
      } else {
        alpha_eq = 0.;
        y_eq = 0.;
        z_eq = 0.;
      }

    }
    else { /* root in [0.5(pmin+pctau),pctau] */
      alpha_eq = alpha2;
      y_eq = y2;
      z_eq = z2;
    }
  }
  else { /* root in [pmin,0.5(pmin+pctau)] */

    alpha_eq = alpha1;
    y_eq = y1;
    z_eq = z1;

  }

  *palpha_eq = alpha_eq;
  *py_eq = y_eq;
  *pz_eq = z_eq;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
