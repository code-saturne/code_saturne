/*============================================================================
 * Libby-Williams gas combustion model.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_physical_constants.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_properties.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"
#include "mesh/cs_mesh.h"
#include "rayt/cs_rad_transfer.h"
#include "turb/cs_turbulence_model.h"

#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_boundary_conditions.h"
#include "cogz/cs_combustion_ht_convert.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_lw.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_lw.cpp
        Libby-Williams gas combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute PDF parameters for 2-point Libby-Williams using the
 *        3rd order moment deduced by recurrences on beta functions.
 *
 * Remark: modified curl hypothesis
 * -------
 *
 * Based on the mean value of a variable, extrema, and the variance of that
 * variable, we deduce 2 states around the mean state and an amplitude for
 * each state.
 *
 * The result is:
 * -------------
 *
 * Compute parameters associated to Dirac functions.
 *
 * The positions of peaks are:
 *   [F[0],Y1[0]] and [F[0],Y1[1]]
 *   [F[1],Y2[0]] and [F[1],Y2[1]]
 * Their respective amplitudes are:
 *   D1[0] and D1[1]
 *   D2[0] and D2[1]
 *
 * \param[in]   ampen1   total amplitude of peaks
 * \param[in]   valmoy   mean value of the variable
 * \param[in]   valvar   variance of the variable
 * \param[in]   valmin   minimum value of the variable
 * \param[in]   valmax   maximum value of the variable
 * \param[out]  exit01   state 1
 * \param[out]  exit02   state 2
 * \param[out]  ampl01   amplitude 1
 * \param[out]  ampl02   amplitude 2
 */
/*----------------------------------------------------------------------------*/

static void
_lwcurl(const cs_real_t   ampen1,
        const cs_real_t   valmoy,
        const cs_real_t   valvar,
        const cs_real_t   valmin,
        const cs_real_t   valmax,
        cs_real_t        &exit01,
        cs_real_t        &exit02,
        cs_real_t        &ampl01,
        cs_real_t        &ampl02)
{
  // Test on total amplitude of peaks. If it is too small
  // we place both on the mean state.

  constexpr cs_real_t epsi = 1.e-6;

  if (ampen1 > epsi) {

    // Test on the variance. If it is too small
    // we place both on the mean state.

    if ((valvar > epsi)) {
      // Adimensionalize for this computation.

      cs_real_t moyadm = (valmoy-valmin) / (valmax-valmin);
      cs_real_t varadm = valvar / cs_math_pow2(valmax-valmin);

      // Adimensional 3rd order moment
      cs_real_t tvvadm = 2. * cs_math_pow2(varadm)
                            * ((1.-2.*moyadm)/((moyadm*(1-moyadm))+varadm));

      // Dimensional 3rd order moment
      cs_real_t tvv = cs_math_pow3(valmax-valmin) * tvvadm;

      // Compute polynomial term for moment 3
      cs_real_t c = 4. + cs_math_pow2(tvv) / cs_math_pow3(valvar);

      // Sign of root
      cs_real_t rsgn = ((1.-moyadm) > moyadm) ? 1. : -1.;
      cs_real_t d = 0.5 + rsgn*sqrt((c-4.)/(4.*c));

      // Amplitudes of peaks
      ampl01 = ampen1 * d;
      ampl02 = ampen1 - ampl01;

      // Positions of peaks
      exit01 = valmoy - sqrt((1.-d)/(d)*valvar);
      exit02 = valmoy + sqrt((d)/(1.-d)*valvar);

      exit01 = std::max(valmin, std::min(exit01, valmax));
      exit02 = std::max(valmin, std::min(exit02, valmax));
    }
    else {
      // Weak variance

      ampl01 = ampen1 * 0.5;
      ampl02 = ampl01;

      exit01 = valmoy;
      exit02 = valmoy;
    }

  }
  else {
    // Weak total amplitude

    ampl01 = ampen1 * 0.5;
    ampl02 = ampl01;

    exit01 = valmoy;
    exit02 = valmoy;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute values of function G.
 *
 * \param[in]   f     mixture fraction
 * \param[in]   fm    mixture fraction mean
 * \param[in]   fm2m  mixture fraction variance
 * \param[in]   yp2m  massfraction variance
 *
 * \return value of g at point f
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_lwcgfu(cs_real_t  f,
        cs_real_t  fm,
        cs_real_t  yfp2m,
        cs_real_t  fp2m)
{
  constexpr cs_real_t epsi = 1.e-9;
  cs_real_t gfunc = 0.;

  if (fp2m <= epsi)
    gfunc = 1.;
  else
    gfunc = (f-fm) * sqrt(1. + yfp2m/fp2m);

  return gfunc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute PDF parameters for 3-point Libby-Williams with curl
 *        or modified curl hypothesis.
 *
 * Remarks:
 * --------
 *
 * In a (F, Yf) diagram, we build 2 lines:
 * - The complete combustion line
 * - The mixing line
 *
 * In this domain, we build 2 peaks on F which defined a 3rd line on
 * which we define a curvilinear abscissa G.
 *
 * The result is:
 * -------------
 *
 * Compute parameters associated to Dirac functions.
 *
 * The positions of peaks are:
 *   [F[0],Y1[0]] and [F[0],Y1[1]]
 *   [F[1],Y2[0]] and [F[1],Y2[1]]
 * Their respective amplitudes are:
 *   D1[0] and D1[1]
 *   D2[0] and D2[1]
 * For each Dirac, compute:
 *   Temperature Ti(j),
 *   Density RHOi(j),
 *   Chemical source term Wi(j),
 *     i being the positiion on F of the Dirac peak
 *     j being the positiion on Yf of the Dirac peak
 *
 * \param[in]  n_cells  number of associated cells
 * \param[in]  fm       mean of mixture fraction
 * \param[in]  fp2m     variance of the mixture fraction
 * \param[in]  yfm      mean of the mass fraction
 * \param[in]  yfp2m    variance of the mass fraction
 */
/*----------------------------------------------------------------------------*/

static void
_pdflwc(cs_lnum_t         n_cells,
        const cs_real_t  *fm,
        const cs_real_t  *fp2m,
        const cs_real_t  *yfm,
        const cs_real_t  *yfp2m)
{
  // Call counter
  static int n_calls = 0;
  n_calls += 1;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  cs_real_t f[CS_COMBUSTION_GAS_MAX_DIRAC], y[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t d[CS_COMBUSTION_GAS_MAX_DIRAC], g[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t h[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t teml[CS_COMBUSTION_GAS_MAX_DIRAC], maml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t w[CS_COMBUSTION_GAS_MAX_DIRAC], rhol[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t theta[CS_COMBUSTION_GAS_MAX_DIRAC];

  const cs_real_t epsi = 1.e-09;

  cs_real_t gmin = HUGE_VALF;
  cs_real_t gmax = -HUGE_VALF;

  // Get variables and coefficients

  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_ym1 = cm->ym[0]->val;
  cs_real_t *cpro_ym2 = cm->ym[1]->val;
  cs_real_t *cpro_ym3 = cm->ym[2]->val;
  cs_real_t *cpro_tsc = cm->lw.tsc->val;
  cs_real_t *cpro_mam = cm->lw.mam->val;

  cs_real_t *cpro_fmel[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_fmal[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_teml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_tscl[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_rhol[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_maml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_ampl[CS_COMBUSTION_GAS_MAX_DIRAC];

  const int n_dirac = cm->lw.n_dirac;
  for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
    cpro_fmel[dirac_id] = cm->lw.fmel[dirac_id]->val;
    cpro_fmal[dirac_id] = cm->lw.fmal[dirac_id]->val;
    cpro_teml[dirac_id] = cm->lw.teml[dirac_id]->val;
    cpro_tscl[dirac_id] = cm->lw.tscl[dirac_id]->val;
    cpro_rhol[dirac_id] = cm->lw.rhol[dirac_id]->val;
    cpro_maml[dirac_id] = cm->lw.maml[dirac_id]->val;
    cpro_ampl[dirac_id] = cm->lw.ampl[dirac_id]->val;
  }

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
    coefg[i] = 0;

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  const cs_real_t p0 = fp->p0;
  const cs_real_t ro0 = fp->ro0;

  const cs_real_t srrom = cm->srrom;
  const cs_real_t ta = cm->lw.ta;
  const cs_real_t tstar = cm->lw.tstar;
  const cs_real_t h_max = cm->lw.hmax;
  const cs_real_t h_min = cm->lw.hmin;
  const cs_real_t coeff1 = cm->lw.coeff1;
  const cs_real_t coeff2 = cm->lw.coeff2;
  const cs_real_t coeff3 = cm->lw.coeff3;
  const cs_real_t vref = cm->lw.vref;
  const cs_real_t lref = cm->lw.lref;
  const cs_real_t *fs = cm->fs;
  const cs_real_t *wmolg = cm->wmolg;

  cs_real_t f_max = cm->lw.fmax;
  cs_real_t f_min = cm->lw.fmin;

  /* Loop on cells
     ------------- */

  const int n_gas_species = cm->n_gas_species;

  bool update_rho = false;
  if (n_calls > 1 || cs_restart_get_field_read_status(CS_F_(rho)->id) == 1)
    update_rho = true;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t fmp    =  std::max(std::min(f_max, fm[c_id]), f_min);
    cs_real_t yfmp   =  std::max(std::min(yfm[c_id], fmp),
                                 (fmp-fs[0]) / (1.-fs[0]));
    yfmp = std::max(yfmp, 0.);
    cs_real_t fp2mp  =  std::max(std::min(fp2m[c_id],
                                          (std::min(f_max,
                                                    (1-fs[0])*yfmp+fs[0])-fmp)
                                          * (fmp-std::max(f_min, yfmp))), 0.);
    cs_real_t yfp2mp
      =  std::max(std::min(yfp2m[c_id],
                           (fmp-yfmp) * (yfmp-std::max(0.,
                                                       (fmp-fs[0])/(1.-fs[0])))),
                  0.);

    // no PDF case

    if (fp2mp < epsi && yfp2mp < epsi) {

      cs_real_t sum1 = 0., sum2 = 0., sum3 = 0., sum4 = 0.;
      cs_real_t sum5 = 0., sum6 = 0., sum15 = 0., sum16 = 0.;

      // For each Dirac peak:
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {

        // Compute f, Y, Amplitude of each DIRAC == mean Val

        d[dirac_id] = 1. / n_dirac;
        f[dirac_id] = fmp;
        y[dirac_id] = yfmp;

        // Compute Enthalpy

        h[dirac_id] =  ((h_max-h_min)*f[dirac_id] + h_min*f_max - h_max*f_min)
                      / (f_max-f_min);

        // Compute mass fraction of gasses (F, O and P)

        cs_real_t yfuel = y[dirac_id];
        cs_real_t yoxyd = 1. - (coeff3+1.)*f[dirac_id] + coeff3*y[dirac_id];
        cs_real_t yprod = 1. - yfuel - yoxyd;
        cs_real_t yo2   = coeff1 - (coeff1 + coeff2) * f[dirac_id]
                                 + coeff2 * y[dirac_id];

        // Compute molar mass and temperature

        coefg[0] = yfuel;
        coefg[1] = yoxyd;
        coefg[2] = yprod;

        cs_real_t nbmol = 0;
        for (int igg = 0; igg < n_gas_species; igg++)
          nbmol += coefg[igg]/wmolg[igg];
        maml[dirac_id] = 1./nbmol;

        // Temperature for each peak

        teml[dirac_id] = cs_gas_combustion_h_to_t(coefg, h[dirac_id]);

        // Density for each peak

        if (update_rho)
          rhol[dirac_id] =   p0 * maml[dirac_id]
                           / (cs_physical_constants_r*teml[dirac_id]);
        else
          rhol[dirac_id] = ro0;

        // Compute source term for scalar YFM for each peak

        theta[dirac_id] =   ta / teml[dirac_id]
                          * (1. - teml[dirac_id] / tstar);
        w[dirac_id] = vref / lref * (- d[dirac_id]*rhol[dirac_id]
                                     * yfuel*yo2
                                     * exp(-theta[dirac_id]));

        // Molar mass of mixture

        sum1 += d[dirac_id]*maml[dirac_id];

        // Mixture temperature

        sum2 += d[dirac_id]*teml[dirac_id];

        // Temperature / Molar mass

        sum3 += d[dirac_id]*teml[dirac_id]/maml[dirac_id];

        // Mass fractions of global species

        sum4 += yfuel*d[dirac_id];
        sum5 += yoxyd*d[dirac_id];
        sum6 += yprod*d[dirac_id];
        sum15 += rhol[dirac_id]*d[dirac_id];
        sum16 += w[dirac_id];

        // Store properties

        cpro_ampl[dirac_id][c_id] = d[dirac_id];
        cpro_fmel[dirac_id][c_id] = f[dirac_id];
        cpro_fmal[dirac_id][c_id] = y[dirac_id];
        cpro_teml[dirac_id][c_id] = teml[dirac_id];
        cpro_maml[dirac_id][c_id] = maml[dirac_id];
        cpro_rhol[dirac_id][c_id] = rhol[dirac_id];
        cpro_tscl[dirac_id][c_id] = w[dirac_id];

      } // Loop on Diracs

      cpro_mam[c_id]  = sum1;
      cpro_temp[c_id] = sum2;
      cpro_ym1[c_id]  = sum4;
      cpro_ym2[c_id]  = sum5;
      cpro_ym3[c_id]  = sum6;
      cpro_tsc[c_id]  = sum16;

      // Mixture density

      if (update_rho) {
        cs_real_t temsmm = sum3;
        crom[c_id] =    srrom * crom[c_id]
                     + (1.-srrom) * (p0/(cs_physical_constants_r*temsmm));
      }

    }

    // PDF case

    else {

      if  (   (fp2mp > epsi && yfp2mp < epsi)
           || (fp2mp < epsi && yfp2mp > epsi)) {

        // case 1: vertical straight line
        if (fp2mp < epsi) {

          // Extrema choices
          gmin = std::max(-yfmp, -(yfmp-(fmp-fs[0]) / (1.-fs[0])));
          gmax = fmp - yfmp;

          // Diracs amplitudes computation
          d[0] = gmax / (gmax-gmin);
          d[1] = (1. - d[0]);

          // Compute variance of curvilinear abscissa (GP2M)
          cs_real_t gp2m = fp2mp + yfp2mp;

          // Test on gp2m
          g[0] = -sqrt(-gmin/gmax*gp2m);
          g[1] = -d[0] * g[0] / d[1];

          for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
            f[dirac_id] = fmp + g[dirac_id] * sqrt(fp2mp/gp2m);
            y[dirac_id] = yfmp + g[dirac_id] * sqrt(yfp2mp/gp2m);
          }

        }

        // case 2 : horizontal straight line

        else if (yfp2mp < epsi) {

          // Computation of the different intersections of the lines
          cs_real_t f12 = yfmp;
          cs_real_t f13 = fs[0] + yfmp * (1-fs[0]);
          cs_real_t f14 = f_min;
          cs_real_t f15 = f_max;

          // Curvilinear coordinate computation
          // call to G function
          cs_real_t g2max = _lwcgfu(f15, fmp, yfp2mp, fp2mp);
          cs_real_t g3max = _lwcgfu(f13, fmp, yfp2mp, fp2mp);
          cs_real_t g1min = _lwcgfu(f12, fmp, yfp2mp, fp2mp);
          cs_real_t g2min = _lwcgfu(f14, fmp, yfp2mp, fp2mp);

          // Extrema choices
          gmin = std::max(g1min, g2min);
          gmax = std::min(g2max, g3max);

          // Computation of Dirac amplitudes
          d[0] = gmax / (gmax-gmin);
          d[1] = (1. - d[0]);

          // Computation of curviline coordinate variance (gp2m)
          cs_real_t gp2m = fp2mp + yfp2mp;

          g[0] = -sqrt(-gmin/gmax*gp2m);
          g[1] = -d[0] * g[0] / d[1];

          for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
            f[dirac_id] = fmp + g[dirac_id] * sqrt(fp2mp/gp2m);
            y[dirac_id] = yfmp + g[dirac_id] * sqrt(yfp2mp/gp2m);
          }

        }

      }

      else if ((fp2mp > epsi) && (yfp2mp > epsi)) {

        // case 3: parallelism and mix line

        if (   ((yfp2mp / fp2mp)  <  1. + epsi)
            && ((yfp2mp / fp2mp)  >  1. - epsi)) {

          cs_real_t aa1 = 1.;
          cs_real_t bb1 = yfmp - fmp;
          cs_real_t aa3 = 1. / (1.-fs[0]);
          cs_real_t bb3 = -fs[0] / (1.-fs[0]);
          cs_real_t aa6 = 0.;
          cs_real_t bb6 = 0.;

          // computation of the different intersections of the lines.
          cs_real_t f13 = (bb3-bb1) / (aa1-aa3);
          cs_real_t f14 = f_min;
          cs_real_t f15 = f_max;
          cs_real_t f16 = (bb6-bb1) / (aa1-aa6);

          // Curvilinear coordinate computation
          // call to G function
          cs_real_t g2max = _lwcgfu(f15, fmp, yfp2mp, fp2mp);
          cs_real_t g3max = _lwcgfu(f13, fmp, yfp2mp, fp2mp);
          cs_real_t g2min = _lwcgfu(f14, fmp, yfp2mp, fp2mp);
          cs_real_t g3min = _lwcgfu(f16, fmp, yfp2mp, fp2mp);

          gmin = std::max(g2min, g3min);
          gmax = std::min(g2max, g3max);

          // computation of Dirac amplitudes
          d[0] = gmax / (gmax-gmin);
          d[1] = (1. - d[0]);

          // Computation of curvilinear coordinate variance (gp2m)
          cs_real_t gp2m = fp2mp + yfp2mp;

          // Test on gp2m
          g[0] = -sqrt(-gmin/gmax*gp2m);
          g[1] =  -d[0] * g[0] / d[1];

          for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
            f[dirac_id] = fmp + g[dirac_id] * sqrt(fp2mp/gp2m);
            y[dirac_id] = yfmp + g[dirac_id] * sqrt(yfp2mp/gp2m);
          }

        }

        // case 4 : parallelism with complete combustion line

        else if (   ((sqrt(yfp2mp/fp2mp) * (1.-fs[0])) < 1. + epsi)
                 && ((sqrt(yfp2mp/fp2mp) * (1.-fs[0])) > 1. - epsi)) {

          cs_real_t aa1 = sqrt( yfp2mp/fp2mp);
          cs_real_t bb1 = yfmp - sqrt( yfp2mp/fp2mp ) * fmp;
          cs_real_t aa2 = 1.;
          cs_real_t bb2 = 0.;
          cs_real_t aa6 = 0.;
          cs_real_t bb6 = 0.;

          // Computation of the different intersections of the lines
          cs_real_t f12 = (bb2-bb1) / (aa1-aa2);
          cs_real_t f14 = f_min;
          cs_real_t f15 = f_max;
          cs_real_t f16 = (bb6-bb1) / (aa1-aa6);

          // Curvilinear coordinate extrema computation
          // call to G function
          cs_real_t g2max = _lwcgfu(f15, fmp, yfp2mp, fp2mp);
          cs_real_t g1max = _lwcgfu(f12, fmp, yfp2mp, fp2mp);
          cs_real_t g2min = _lwcgfu(f14, fmp, yfp2mp, fp2mp);
          cs_real_t g3min = _lwcgfu(f16, fmp, yfp2mp, fp2mp);

          gmin = std::max(g2min, g3min);
          gmax = std::min(g1max, g2max);

          // Dirac amplitudes computations
          d[0] = gmax / (gmax-gmin);
          d[1] = (1. - d[0]);

          // ---> computation of curviline coordinate variance (gp2m)

          cs_real_t gp2m = fp2mp + yfp2mp;

          // Test on gp2m

          g[0] = -sqrt(-gmin/gmax*gp2m);
          g[1] = -d[0] * g[0] / d[1];

          for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
            f[dirac_id] = fmp + g[dirac_id] * sqrt(fp2mp/gp2m);
            y[dirac_id] = yfmp + g[dirac_id] * sqrt(yfp2mp/gp2m);
          }

        }

        // General case

        else {

          // Expression for the different lines: Y=AX+B

          cs_real_t aa1 = -sqrt(yfp2mp/fp2mp);
          cs_real_t bb1 = yfmp - aa1 * fmp;
          cs_real_t aa2 = 1.;
          cs_real_t bb2 = 0.;
          cs_real_t aa3 = 1. / (1.-fs[0]);
          cs_real_t bb3 = -fs[0] / (1.-fs[0]);
          cs_real_t aa6 = 0.;
          cs_real_t bb6 = 0.;

          // Computation of the different intersections of the lines.

          cs_real_t f12 = (bb2-bb1) / (aa1-aa2);
          cs_real_t f13 = (bb3-bb1) / (aa1-aa3);
          cs_real_t f14 = f_min;
          cs_real_t f15 = f_max;
          cs_real_t f16 = (bb6-bb1) / (aa1-aa6);

          // Curvilinear coordinate extrema computation.
          // call to G function

          cs_real_t g1max = _lwcgfu(f12, fmp, yfp2mp, fp2mp);
          cs_real_t g2max = _lwcgfu(f15, fmp, yfp2mp, fp2mp);
          cs_real_t g3max = _lwcgfu(f13, fmp, yfp2mp, fp2mp);
          cs_real_t g1min = _lwcgfu(f12, fmp, yfp2mp, fp2mp);
          cs_real_t g2min = _lwcgfu(f14, fmp, yfp2mp, fp2mp);
          cs_real_t g3min = _lwcgfu(f16, fmp, yfp2mp, fp2mp);
          cs_real_t g4min = _lwcgfu(f13, fmp, yfp2mp, fp2mp);

          /* Computation of the parameters of the two Dirac peaks
             ---------------------------------------------------- */

          // Extrema choice according to the lines slope
          if (aa1 > aa3) {
            gmin = std::max (std::max(g2min, g3min), g4min);
            gmax = std::min (g1max, g2max);
          }
          else if ((aa1 > aa2) && (aa1 <= aa3)) {
            gmin = std::max(g2min, g3min);
            gmax = std::min(std::min(g1max, g2max), g3max);
          }
          else if (aa1 <= aa2) {
            gmin = std::max(std::max(g1min, g2min), g3min);
            gmax = std::min(g2max, g3max);
          }

          // Dirac amplitudes computations
          d[0] = gmax / (gmax-gmin);
          d[1] = (1. - d[0]);

          // Compute variance of curvilinear abscissa (GP2M)
          cs_real_t gp2m = fp2mp + yfp2mp;

          g[0] = -sqrt(gmin/gmax*gp2m);
          g[1] = -d[0] * g[0] / d[1];

          for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
            f[dirac_id] = fmp + g[dirac_id] * sqrt( fp2mp/gp2m);
            y[dirac_id] = yfmp + g[dirac_id] * sqrt( yfp2mp/gp2m);
          }

        }
      }

      /* Computations of the thermochemical quantities at the two peaks
         -------------------------------------------------------------- */

      // Enthalpy computation in 1 and 2
      cs_real_t sum7  = 0., sum8  = 0., sum9  = 0., sum10 = 0.;
      cs_real_t sum11 = 0., sum12 = 0., sum17 = 0.;

      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {

        h[dirac_id] =   ((h_max-h_min)*f[dirac_id] + h_min*f_max - h_max*f_min)
                      / (f_max-f_min);

        // Computation of the mass fraction of the gas (F, O et P) in 1 and 2.
        cs_real_t yfuel = y[dirac_id];
        cs_real_t yoxyd = 1. - (coeff3+1.)*f[dirac_id] + coeff3*y[dirac_id];
        cs_real_t yprod = 1. - yfuel - yoxyd;
        cs_real_t yo2   = coeff1 - (coeff1 + coeff2) * f[dirac_id]
                                  + coeff2 * y[dirac_id];

        // Computation of molar mass and temperature in 1 and 2.
        coefg[0] = yfuel;
        coefg[1] = yoxyd;
        coefg[2] = yprod;

        // Molar mass for peaks 1 and 2
        cs_real_t nbmol = 0.;
        for (int igg = 0; igg < n_gas_species; igg++)
          nbmol += coefg[igg]/wmolg[igg];
        maml[dirac_id] = 1./nbmol;

        // Computation of temperature for peaks 1 and 2.
        teml[dirac_id] = cs_gas_combustion_h_to_t(coefg, h[dirac_id]);

        // Computation of density in 1 and 2
        if (update_rho)
          rhol[dirac_id] =   p0 * maml[dirac_id]
                           / (cs_physical_constants_r*teml[dirac_id]);
        else
          rhol[dirac_id] = ro0;

        // Source term computation in 1 and 2 for scalar yfm

        theta[dirac_id] =   ta / teml[dirac_id]
                          * (1. - teml[dirac_id] / tstar);

        w[dirac_id] = vref / lref * (  -d[dirac_id]*rhol[dirac_id]
                                     * yfuel*yo2
                                     * exp( -theta[dirac_id]));

        // Mix molar mass
        sum7 += d[dirac_id]*maml[dirac_id];

        // Mix temperature
        sum8 += d[dirac_id]*teml[dirac_id];

        // Temperature / molar mass
        sum9 += d[dirac_id]*teml[dirac_id]/maml[dirac_id];

        // Global species mass fractions
        sum10 += yfuel*d[dirac_id];
        sum11 += yoxyd*d[dirac_id];
        sum12 += yprod*d[dirac_id];
        sum17 += w[dirac_id];

        // Property storage in fields
        cpro_ampl[dirac_id][c_id] = d[dirac_id];
        cpro_fmel[dirac_id][c_id] = f[dirac_id];
        cpro_fmal[dirac_id][c_id] = y[dirac_id];
        cpro_maml[dirac_id][c_id] = maml[dirac_id];
        cpro_teml[dirac_id][c_id] = teml[dirac_id];
        cpro_rhol[dirac_id][c_id] = rhol[dirac_id];
        cpro_tscl[dirac_id][c_id] = w[dirac_id];

      }

      cpro_mam[c_id]  = sum7;
      cpro_temp[c_id] = sum8;
      cpro_ym1[c_id]  = sum10;
      cpro_ym2[c_id]  = sum11;
      cpro_ym3[c_id]  = sum12;
      cpro_tsc[c_id]  = sum17;

      if (update_rho) {
        cs_real_t temsmm = sum9;
        crom[c_id] =    srrom * crom[c_id]
                     + (1.-srrom) * (p0/(cs_physical_constants_r*temsmm));
      }

    }

  } // End of loop on cells
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute PDF parameters for 3-point Libby-Williams with curl
 *        or modified curl hypothesis.
 *
 * Remarks:
 * --------
 *
 * In a (F, Yf) diagram, we build 2 lines:
 * - The complete combustion line
 * - The mixing line
 *
 * In this domain, we build 2 peaks on F which defined a 3rd line on
 * which we define a curvilinear abscissa G.
 *
 * The result is:
 * -------------
 *
 * Compute parameters associated to Dirac functions.
 *
 * The positions of peaks are:
 *   [F[0],Y1[0]] and [F[0],Y1[1]]
 *   [F[1],Y2[0]] and [F[1],Y2[1]]
 * Their respective amplitudes are:
 *   D1[0] and D1[1]
 *   D2[0] and D2[1]
 * For each Dirac, compute:
 *   Temperature Ti(j),
 *   Density RHOi(j),
 *   Chemical source term Wi(j),
 *     i being the positiion on F of the Dirac peak
 *     j being the positiion on Yf of the Dirac peak
 *
 * \param[in]  n_cells  number of associated cells
 * \param[in]  fm       mean of mixture fraction
 * \param[in]  fp2m     variance of the mixture fraction
 * \param[in]  yfm      mean of the mass fraction
 * \param[in]  yfp2m    variance of the mass fraction
 * \param[in]  coyfp    covariance
 */
/*----------------------------------------------------------------------------*/

static void
_pdfpp3(const cs_lnum_t   n_cells,
        const cs_real_t  *fm,
        const cs_real_t  *fp2m,
        const cs_real_t  *yfm,
        const cs_real_t  *yfp2m,
        const cs_real_t  *coyfp)
{
  // Call counter
  static int n_calls = 0;
  n_calls += 1;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  cs_real_t f[CS_COMBUSTION_GAS_MAX_DIRAC], y[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t d[CS_COMBUSTION_GAS_MAX_DIRAC], h[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t teml[CS_COMBUSTION_GAS_MAX_DIRAC], maml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t w[CS_COMBUSTION_GAS_MAX_DIRAC], rhol[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t theta[CS_COMBUSTION_GAS_MAX_DIRAC];

  const cs_real_t epsi = 1.e-6;

  // Get variables and coefficients

  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_ym1 = cm->ym[0]->val;
  cs_real_t *cpro_ym2 = cm->ym[1]->val;
  cs_real_t *cpro_ym3 = cm->ym[2]->val;
  cs_real_t *cpro_tsc = cm->lw.tsc->val;
  cs_real_t *cpro_mam = cm->lw.mam->val;

  cs_real_t *cpro_fmel[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_fmal[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_teml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_tscl[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_rhol[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_maml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_ampl[CS_COMBUSTION_GAS_MAX_DIRAC];

  const int n_dirac = cm->lw.n_dirac;
  for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
    cpro_fmel[dirac_id] = cm->lw.fmel[dirac_id]->val;
    cpro_fmal[dirac_id] = cm->lw.fmal[dirac_id]->val;
    cpro_teml[dirac_id] = cm->lw.teml[dirac_id]->val;
    cpro_tscl[dirac_id] = cm->lw.tscl[dirac_id]->val;
    cpro_rhol[dirac_id] = cm->lw.rhol[dirac_id]->val;
    cpro_maml[dirac_id] = cm->lw.maml[dirac_id]->val;
    cpro_ampl[dirac_id] = cm->lw.ampl[dirac_id]->val;
  }

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
    coefg[i] = 0;

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  const cs_real_t p0 = fp->p0;
  const cs_real_t ro0 = fp->ro0;

  const cs_real_t srrom = cm->srrom;
  const cs_real_t ta = cm->lw.ta;
  const cs_real_t tstar = cm->lw.tstar;
  const cs_real_t h_max = cm->lw.hmax;
  const cs_real_t h_min = cm->lw.hmin;
  const cs_real_t coeff1 = cm->lw.coeff1;
  const cs_real_t coeff2 = cm->lw.coeff2;
  const cs_real_t coeff3 = cm->lw.coeff3;
  const cs_real_t vref = cm->lw.vref;
  const cs_real_t lref = cm->lw.lref;
  const cs_real_t *fs = cm->fs;
  const cs_real_t *wmolg = cm->wmolg;

  cs_real_t f_max = cm->lw.fmax;
  cs_real_t f_min = cm->lw.fmin;

  // Counters

  cs_gnum_t clif = 0, cliy = 0, clifp2 = 0, clyfp2 = 0, clicoy = 0;
  cs_gnum_t cliy1 = 0, cliy2 = 0, cliy2p = 0;

  cs_real_t ymin[2], ymax[2], y2p[2];

  /* Loop on cells
     ------------- */

  const int n_gas_species = cm->n_gas_species;

  bool update_rho = false;
  if (n_calls > 1 || cs_restart_get_field_read_status(CS_F_(rho)->id) == 1)
    update_rho = true;
  else {
    f_min = 4.405286343612334e-02;
    f_max = 5.506607929515418e-02;
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    // F
    cs_real_t fmp = fm[c_id];
    if ((fmp <= f_min) || (fmp >= f_max)) {
      fmp = std::max(std::min(f_max, fmp), f_min);
      clif += 1;
    }

    // Y
    cs_real_t yfmp = yfm[c_id];
    cs_real_t climax = (f_max-fmp)*f_min / (f_max-f_min);
    cs_real_t climin = std::max(0., (fmp-fs[0])/(1.-fs[0]));
    if ((yfmp >= climax) || (yfmp < climin)) {
      yfmp = std::max(climin, std::min(yfmp, climax));
      cliy += 1;
    }

    // FP2M
    cs_real_t fp2mp = fp2m[c_id];
    climax = (f_max -fmp)*(fmp-f_min);
    climin = 0.;
    if ((fp2mp >= climax) || (fp2mp < climin)) {
      fp2mp = std::max(climin, std::min(fp2m[c_id], climax));
      clifp2 += 1;
    }

    // YFP2M
    // yf_max = f_min in the Moreau case
    cs_real_t yfp2mp = yfp2m[c_id];
    climax = (f_min-yfmp)*yfmp;
    climin = 0.;
    if ((yfp2mp >= climax) || (yfp2mp < climin)) {
      yfp2mp = std::max(climin, std::min(yfp2mp, climax));
      clyfp2 += 1;
    }

    // Clip for covariance

    cs_real_t coyfpp = coyfp[c_id];
    climax = sqrt(fp2mp*yfp2mp);
    climin = -sqrt(fp2mp*yfp2mp);
    if (coyfpp >= climax) {
      coyfpp = climax;
      clicoy += 1;
    }
    else if (coyfpp <= climin) {
      coyfpp = climin;
      clicoy += 1;
    }

    cs_real_t yfmpmx = (f_max - fmp)*f_min/(f_max - f_min);
    cs_real_t dyfmp  = (yfmpmx - yfmp);

    /* Test whether we go through the PDF
       ----------------------------------
       (do not go through PDF if no variance either in f or y:
       we could go through the pdf in all cases, as zero variances is handled */

    // no PDF case

    if (   (fp2mp < epsi && yfp2mp < epsi)
        || (yfmp < epsi || dyfmp < epsi)
        || (fmp -f_min < epsi)
        || (fp2mp < epsi)
        || (f_max -fmp < epsi)) {

      cs_real_t sum1 = 0., sum2 = 0., sum3 = 0., sum4 = 0.;
      cs_real_t sum5 = 0., sum6 = 0., sum15 = 0., sum16 = 0.;

      // For each Dirac peak:
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {

        // Compute f, Y, Amplitude of each DIRAC == mean Val

        d[dirac_id] = 1. / n_dirac;
        f[dirac_id] = fmp;
        y[dirac_id] = yfmp;

        // Compute Enthalpy

        h[dirac_id] =  ((h_max-h_min)*f[dirac_id] + h_min*f_max - h_max*f_min)
                      / (f_max-f_min);

        // Compute mass fraction of gasses (F, O and P)

        cs_real_t yfuel = y[dirac_id];
        cs_real_t yoxyd = 1. - (coeff3+1.)*f[dirac_id] + coeff3*y[dirac_id];
        cs_real_t yprod = 1. - yfuel - yoxyd;
        cs_real_t yo2   = coeff1 - (coeff1 + coeff2) * f[dirac_id]
                                 + coeff2 * y[dirac_id];

        // Compute molar mass and temperature

        coefg[0] = yfuel;
        coefg[1] = yoxyd;
        coefg[2] = yprod;

        cs_real_t nbmol = 0;
        for (int igg = 0; igg < n_gas_species; igg++)
          nbmol += coefg[igg]/wmolg[igg];
        maml[dirac_id] = 1./nbmol;

        // Temperature for each peak

        teml[dirac_id] = cs_gas_combustion_h_to_t(coefg, h[dirac_id]);

        // Density for each peak

        if (update_rho)
          rhol[dirac_id] =   p0 * maml[dirac_id]
                           / (cs_physical_constants_r*teml[dirac_id]);
        else
          rhol[dirac_id] = ro0;

        // Compute source term for scalar YFM for each peak

        theta[dirac_id] =   ta / teml[dirac_id]
                          * (1. - teml[dirac_id] / tstar);
        w[dirac_id] = vref / lref * (- d[dirac_id]*rhol[dirac_id]
                                     * yfuel*yo2
                                     * exp(-theta[dirac_id]));

        // BO 27/06 Control sign of W
        // FIXME: check this

        w[dirac_id] = std::min(w[dirac_id], 0.);

        // Molar mass of mixture

        sum1 += d[dirac_id]*maml[dirac_id];

        // Mixture temperature

        sum2 += d[dirac_id]*teml[dirac_id];

        // Temperature / Molar mass

        sum3 += d[dirac_id]*teml[dirac_id]/maml[dirac_id];

        // Mass fractions of global species

        sum4 += yfuel*d[dirac_id];
        sum5 += yoxyd*d[dirac_id];
        sum6 += yprod*d[dirac_id];
        sum15 += rhol[dirac_id]*d[dirac_id];
        sum16 += w[dirac_id];

        // Store properties

        cpro_ampl[dirac_id][c_id] = d[dirac_id];
        cpro_fmel[dirac_id][c_id] = f[dirac_id];
        cpro_fmal[dirac_id][c_id] = y[dirac_id];
        cpro_teml[dirac_id][c_id] = teml[dirac_id];
        cpro_maml[dirac_id][c_id] = maml[dirac_id];
        cpro_rhol[dirac_id][c_id] = rhol[dirac_id];
        cpro_tscl[dirac_id][c_id] = w[dirac_id];

      } // Loop on Diracs

      cpro_mam[c_id]  = sum1;
      cpro_temp[c_id] = sum2;
      cpro_ym1[c_id]  = sum4;
      cpro_ym2[c_id]  = sum5;
      cpro_ym3[c_id]  = sum6;
      cpro_tsc[c_id]  = sum16;

      // Mixture density

      if (update_rho) {
        cs_real_t temsmm = sum3;
        crom[c_id] =    srrom * crom[c_id]
                     + (1.-srrom) * (p0/(cs_physical_constants_r*temsmm));
      }

    }

    // PDF case

    else {

      // Compute F1 and F2 with Curl at F

      cs_real_t f1, f2, cstfa1, cstfa2;

      _lwcurl(1.0, fmp, fp2mp, f_min, f_max, f1, f2, cstfa1, cstfa2);

      // Compute conditional means Y1, Y2

      cs_real_t y2 = ((fmp*yfmp + coyfpp) - f1*yfmp) / (cstfa1*(f2 - f1));
      cs_real_t y1 = (yfmp - cstfa2*y2)/cstfa1;

      ymin[0] = std::max(0., ((f1- fs[0])/(1.-fs[0])));
      ymax[0] = (f_max - f1)*f_min/(f_max - f_min);
      ymin[1] = std::max(0., ((f2- fs[0])/(1.-fs[0])));
      ymax[1] = (f_max - f2)*f_min/(f_max - f_min);

      // clipping for conditional means

      if (y1 >= ymax[0]) {
        y1 = ymax[0];
        cliy1 += 1;
      }
      else if (y1 <= ymin[0]) {
        y1 = ymin[0];
        cliy1 += 1;
      }

      if (y2 >= ymax[1]) {
        y2 = ymax[1];
        cliy2 += 1;
      }
      else if (y2 <= ymin[1]) {
        y2 = ymin[1];
        cliy2 += 1;
      }

      y2p[0]  =   ((cs_math_pow2(yfmp) + yfp2mp) - cstfa2*(cs_math_pow2(y2)))
                / cstfa1
                - cs_math_pow2(y1);

      // clipping for conditional variances

      climax = ((y1-ymin[0])*(ymax[0] - y1));
      climin = 0.;
      if (y2p[0] >= climax) {
        y2p[0] = climax;
        cliy2p += 1;
      }
      else if (y2p[0] <= climin) {
        y2p[0] = climin;
        cliy2p += 1;
      }

      _lwcurl(cstfa1, y1, y2p[0], ymin[0], ymax[0],
              y[0], y[1], d[0], d[1]);

      // Dirac parameters at F1

      f[0] = f1;
      f[1] = f1;

      // Dirac parameters at F2

      f[2] = f2;
      y[2] = y2;
      d[2] = cstfa2;

      // Compute enthalpies at 1 and 2

      cs_real_t sum1  = 0., sum2  = 0., sum3  = 0., sum4 = 0.;
      cs_real_t sum5 = 0., sum6 = 0., sum16 = 0.;

      // Loop on each Dirac.

      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {

        // Compute Enthalpye

        h[dirac_id] = (  (h_max-h_min) * f[dirac_id]
                       + h_min*f_max - h_max*f_min) / (f_max-f_min);

        // Compute mass fraction of gases (F, O and P) for each peak.

        cs_real_t yfuel = y[dirac_id];
        cs_real_t yoxyd = 1. - (coeff3+1.)*f[dirac_id] + coeff3*y[dirac_id];
        cs_real_t yprod = 1. - yfuel - yoxyd;
        cs_real_t yo2   =   coeff1 - (coeff1 + coeff2) * f[dirac_id]
                          + coeff2 * y[dirac_id];

        // Compute molar mass and temperature for each peak.

        coefg[0] = yfuel;
        coefg[1] = yoxyd;
        coefg[2] = yprod;

        // Molar mass for each peak.

        cs_real_t nbmol = 0;
        for (int igg = 0; igg < n_gas_species; igg++)
          nbmol += coefg[igg]/wmolg[igg];
        maml[dirac_id] = 1./nbmol;

        // Temperature for each peak

        teml[dirac_id] = cs_gas_combustion_h_to_t(coefg, h[dirac_id]);

        // Density for each peak

        if (update_rho)
          rhol[dirac_id] =   p0 * maml[dirac_id]
                           / (cs_physical_constants_r*teml[dirac_id]);
        else
          rhol[dirac_id] = ro0;

        // Compute source term for scalar YFM for each peak.

        theta[dirac_id] =   ta / teml[dirac_id]
                          * (1. - teml[dirac_id] / tstar);

        w[dirac_id] = vref / lref * (- d[dirac_id]*rhol[dirac_id]
                                     * yfuel*yo2
                                     * exp(-theta[dirac_id]));

        // Control sign of W

        w[dirac_id] = std::min(w[dirac_id], 0.);

        // Molar mass of mixture

        sum1 += d[dirac_id]*maml[dirac_id];

        // Mixture temperature

        sum2 += d[dirac_id]*teml[dirac_id];

        // Temperature / Molar mass

        sum3 += d[dirac_id]*teml[dirac_id]/maml[dirac_id];

        // Mass fractions of global species

        sum4 += yfuel*d[dirac_id];
        sum5 += yoxyd*d[dirac_id];
        sum6 += yprod*d[dirac_id];
        sum16 += w[dirac_id];

        // Store properties

        cpro_ampl[dirac_id][c_id] = d[dirac_id];
        cpro_fmel[dirac_id][c_id] = f[dirac_id];
        cpro_fmal[dirac_id][c_id] = y[dirac_id];
        cpro_maml[dirac_id][c_id] = maml[dirac_id];
        cpro_teml[dirac_id][c_id] = teml[dirac_id];
        cpro_rhol[dirac_id][c_id] = rhol[dirac_id];
        cpro_tscl[dirac_id][c_id] = w[dirac_id];

      } // Loop on Diracs

      cpro_mam[c_id]  = sum1;
      cpro_temp[c_id] = sum2;
      cpro_ym1[c_id]  = sum4;
      cpro_ym2[c_id]  = sum5;
      cpro_ym3[c_id]  = sum6;
      cpro_tsc[c_id]  = sum16;

      // Mixture density

      if (update_rho) {
        cs_real_t temsmm = sum3;
        crom[c_id] =    srrom * crom[c_id]
                     + (1.-srrom) * (p0/(cs_physical_constants_r*temsmm));
      }

    } // End if PDF case

  } // End of loop on cells
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute PDF parameters for 4-point Libby-Williams with curl
 *        or modified curl hypothesis.
 *
 * Remarks:
 * --------
 *
 * In a (F, Yf) diagram, we build 2 lines:
 * - The complete combustion line
 * - The mixing line
 *
 * In this domain, we build 2 peaks on F which will be each be doubled, with
 * a Curl on Yf.
 *
 * The result is:
 * -------------
 *
 * Compute parameters associated to Dirac functions.
 *
 * The positions of peaks are:
 *   [F[0],Y1[0]] and [F[0],Y1[1]]
 *   [F[1],Y2[0]] and [F[1],Y2[1]]
 * Their respective amplitudes are:
 *   D1[0] and D1[1]
 *   D2[0] and D2[1]
 * For each Dirac, compute:
 *   Temperature Ti(j),
 *   Density RHOi(j),
 *   Chemical source term Wi(j),
 *     i being the positiion on F of the Dirac peak
 *     j being the positiion on Yf of the Dirac peak
 *
 * \param[in]  n_cells  number of associated cells
 * \param[in]  fm       mean of mixture fraction
 * \param[in]  fp2m     variance of the mixture fraction
 * \param[in]  yfm      mean of the mass fraction
 * \param[in]  yfp2m    variance of the mass fraction
 * \param[in]  coyfp    covariance
 */
/*----------------------------------------------------------------------------*/

static void
_pdfpp4(const cs_lnum_t   n_cells,
        const cs_real_t  *fm,
        const cs_real_t  *fp2m,
        const cs_real_t  *yfm,
        const cs_real_t  *yfp2m,
        const cs_real_t  *coyfp)
{
  // Call counter
  static int n_calls = 0;
  n_calls += 1;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  cs_real_t f[CS_COMBUSTION_GAS_MAX_DIRAC], y[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t d[CS_COMBUSTION_GAS_MAX_DIRAC], h[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t teml[CS_COMBUSTION_GAS_MAX_DIRAC], maml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t w[CS_COMBUSTION_GAS_MAX_DIRAC], rhol[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t theta[CS_COMBUSTION_GAS_MAX_DIRAC];

  const cs_real_t epsi = 1.e-10;

  // Get variables and coefficients

  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_ym1 = cm->ym[0]->val;
  cs_real_t *cpro_ym2 = cm->ym[1]->val;
  cs_real_t *cpro_ym3 = cm->ym[2]->val;
  cs_real_t *cpro_tsc = cm->lw.tsc->val;
  cs_real_t *cpro_mam = cm->lw.mam->val;

  cs_real_t *cpro_fmel[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_fmal[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_teml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_tscl[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_rhol[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_maml[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_ampl[CS_COMBUSTION_GAS_MAX_DIRAC];

  const int n_dirac = cm->lw.n_dirac;
  for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
    cpro_fmel[dirac_id] = cm->lw.fmel[dirac_id]->val;
    cpro_fmal[dirac_id] = cm->lw.fmal[dirac_id]->val;
    cpro_teml[dirac_id] = cm->lw.teml[dirac_id]->val;
    cpro_tscl[dirac_id] = cm->lw.tscl[dirac_id]->val;
    cpro_rhol[dirac_id] = cm->lw.rhol[dirac_id]->val;
    cpro_maml[dirac_id] = cm->lw.maml[dirac_id]->val;
    cpro_ampl[dirac_id] = cm->lw.ampl[dirac_id]->val;
  }

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
    coefg[i] = 0;

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  const cs_real_t p0 = fp->p0;
  const cs_real_t ro0 = fp->ro0;

  const cs_real_t srrom = cm->srrom;
  const cs_real_t ta = cm->lw.ta;
  const cs_real_t tstar = cm->lw.tstar;
  const cs_real_t h_max = cm->lw.hmax;
  const cs_real_t h_min = cm->lw.hmin;
  const cs_real_t f_max = cm->lw.fmax;
  const cs_real_t f_min = cm->lw.fmin;
  const cs_real_t coeff1 = cm->lw.coeff1;
  const cs_real_t coeff2 = cm->lw.coeff2;
  const cs_real_t coeff3 = cm->lw.coeff3;
  const cs_real_t vref = cm->lw.vref;
  const cs_real_t lref = cm->lw.lref;
  const cs_real_t *fs = cm->fs;
  const cs_real_t *wmolg = cm->wmolg;

  // Counters

  cs_gnum_t clfp2m = 0, clyf21 = 0, clcyf1 = 0, clcyf2 = 0;
  cs_gnum_t cly2p1 = 0, cly2p2 = 0, cly2p3 = 0, cly2p4 = 0;
  cs_gnum_t icpt1 = 0, icpt2 = 0;

  cs_real_t wmax = -1.e+10, wmin = 1.e+10, tetmin = 1.e+10, tetmax = -1.e+10;
  cs_real_t o2max = -1.e+10, o2min = 1.e+10;

  cs_real_t mxcfp = -1.e+10, mxcyfp = -1.e+10;
  cs_real_t mxccyf = -1.e+10, mnccyf = -1.e+10;

  cs_real_t maxfp2 = -1.e+10, mxyfp2 = -1.e+10;
  cs_real_t mxcoyf = -1.e+10, mncoyf = -1.e+10;
  cs_real_t mcy2p1 = -1.e+10;

  cs_real_t mcy2p3 = -1.e+10, mcy2p2 = -1.e+10, mcy2p4 = -1.e+10;
  cs_real_t my2p1 = -1.e+10, my2p3 = -1.e+10, my2p2 = -1.e+10, my2p4 = -1.e+10;

  cs_real_t ymin[2], ymax[2], y2p[2], y2pmin[2], y2pmax[2];

  /* Loop on cells
     ------------- */

  const int n_gas_species = cm->n_gas_species;

  bool update_rho = false;
  if (n_calls > 1 || cs_restart_get_field_read_status(CS_F_(rho)->id) == 1)
    update_rho = true;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t fmp = fm[c_id];
    cs_real_t yfmp = std::max(yfm[c_id], 0.);
    cs_real_t fp2mp = fp2m[c_id];
    cs_real_t yfp2mp = yfp2m[c_id];
    cs_real_t coyfpp = coyfp[c_id];

    /* Test whether we go through the PDF
       ----------------------------------
       (do not go through PDF if no variance either in f or y:
       we could go through the pdf in all cases, as zero variances is handled */

    // no PDF case

    if ((fp2mp < epsi) && (yfp2mp < epsi)) {

      icpt1++;

      cs_real_t sum1 = 0., sum2 = 0., sum3 = 0., sum4 = 0.;
      cs_real_t sum5 = 0., sum6 = 0., sum15 = 0., sum16 = 0.;

      // For each Dirac peak:
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {

        // Compute f, Y, Amplitude of each DIRAC == mean Val

        d[dirac_id] = 1. / n_dirac;
        f[dirac_id] = fmp;
        y[dirac_id] = yfmp;

        // Compute Enthalpy

        h[dirac_id] =  ((h_max-h_min)*f[dirac_id] + h_min*f_max - h_max*f_min)
                      / (f_max-f_min);

        // Compute mass fraction of gasses (F, O and P)

        cs_real_t yfuel = y[dirac_id];
        cs_real_t yoxyd = 1. - (coeff3+1.)*f[dirac_id] + coeff3*y[dirac_id];
        cs_real_t yprod = 1. - yfuel - yoxyd;
        cs_real_t yo2   = coeff1 - (coeff1 + coeff2) * f[dirac_id]
                                 + coeff2 * y[dirac_id];

        // Compute molar mass and temperature

        coefg[0] = yfuel;
        coefg[1] = yoxyd;
        coefg[2] = yprod;

        cs_real_t nbmol = 0;
        for (int igg = 0; igg < n_gas_species; igg++)
          nbmol += coefg[igg]/wmolg[igg];
        maml[dirac_id] = 1./nbmol;

        // Temperature for each peak
        teml[dirac_id] = cs_gas_combustion_h_to_t(coefg, h[dirac_id]);

        // Density for each peak
        if (update_rho)
          rhol[dirac_id] =   p0 * maml[dirac_id]
                           / (cs_physical_constants_r*teml[dirac_id]);
        else
          rhol[dirac_id] = ro0;

        // Compute source term for scalar YFM for each peak
        theta[dirac_id] =   ta / teml[dirac_id]
                          * (1. - teml[dirac_id] / tstar);
        w[dirac_id] = vref / lref * (- d[dirac_id]
                                     * yfuel*yo2
                                     * exp(-theta[dirac_id]));

        // FIXME: in pdfpp3, d[dirac_id] above is multiplied by
        // rhol[dirac_id], but not here. Is this normal, or a bug ?

        // Control sign of W
        w[dirac_id] = std::min(w[dirac_id], 0.);

        // Molar mass of mixture
        sum1 += d[dirac_id]*maml[dirac_id];

        // Mixture temperature
        sum2 += d[dirac_id]*teml[dirac_id];

        // Temperature / Molar mass

        sum3 += d[dirac_id]*teml[dirac_id]/maml[dirac_id];

        // Mass fractions of global species

        sum4 += yfuel*d[dirac_id];
        sum5 += yoxyd*d[dirac_id];
        sum6 += yprod*d[dirac_id];
        sum15 += rhol[dirac_id]*d[dirac_id];
        sum16 += w[dirac_id];

        // Store properties

        cpro_ampl[dirac_id][c_id] = d[dirac_id];
        cpro_fmel[dirac_id][c_id] = f[dirac_id];
        cpro_fmal[dirac_id][c_id] = y[dirac_id];
        cpro_teml[dirac_id][c_id] = teml[dirac_id];
        cpro_maml[dirac_id][c_id] = maml[dirac_id];
        cpro_rhol[dirac_id][c_id] = rhol[dirac_id];
        cpro_tscl[dirac_id][c_id] = w[dirac_id];

      } // Loop on Diracs

      cpro_mam[c_id]  = sum1;
      cpro_temp[c_id] = sum2;
      cpro_ym1[c_id]  = sum4;
      cpro_ym2[c_id]  = sum5;
      cpro_ym3[c_id]  = sum6;
      cpro_tsc[c_id]  = sum16;

      // Mixture density

      if (update_rho) {
        cs_real_t temsmm = sum3;
        crom[c_id] =    srrom * crom[c_id]
                     + (1.-srrom) * (p0/(cs_physical_constants_r*temsmm));
      }

    }

    // PDF case

    else {

      icpt2++;

      // Clipping on variance in f

      cs_real_t vfmx = (f_max-fmp)*(fmp-f_min);

      if (fp2mp > vfmx) {
        if ((fp2mp-vfmx) > mxcfp) {
          mxcfp = (fp2mp-vfmx);
          maxfp2 = vfmx;
        }
        fp2mp = vfmx;
        clfp2m += 1;
      }

      // Compute positions and amplitudes at F of peaks with lWCURL at F

      cs_real_t cst = 1.;
      cs_real_t f1, f2, cstfa1, cstfa2;

      _lwcurl(cst, fmp, fp2mp, f_min, f_max,
              f1, f2, cstfa1, cstfa2);

      f[0] = f1;
      f[1] = f1;
      f[2] = f2;
      f[3] = f2;

      // Determine max and min of Yf in F1 et F2

      ymin[0] = std::max(0., ((f1- fs[0])/(1.-fs[0])));
      ymax[0] = f1;
      ymin[1] = std::max(0., ((f2- fs[0])/(1.-fs[0])));
      ymax[1] = f2;

      // covariance clipping based on bounds of conditionl means

      cs_real_t vcyfmx = std::min(ymax[1]*cstfa2*(f2-f1)-yfmp*(fmp-f1),
                                  yfmp*(f2-fmp)-ymin[0]*cstfa1*(f2-f1));

      cs_real_t vcyfmn = std::max(ymin[1]*cstfa2*(f2-f1)-yfmp*(fmp-f1),
                                  yfmp*(f2-fmp)-ymax[0]*cstfa1*(f2-f1));

      if (coyfpp > vcyfmx) {
        if ((coyfpp-vcyfmx) > mxccyf) {
          mxccyf = (coyfpp-vcyfmx);
          mxcoyf = vcyfmx;
        }
        coyfpp = vcyfmx;
        clcyf1 += 1;
      }
      else if (coyfpp < vcyfmn) {
        if ((vcyfmn-coyfpp) > mnccyf) {
          mnccyf = (vcyfmn-coyfpp);
          mncoyf = vcyfmn;
        }
        coyfpp = vcyfmn;
        clcyf2 += 1;
      }

      // Compute conditionnal means Y1, Y2

      cs_real_t y1, y2;
      if ((f2-f1) > epsi) {
        y2 =   (yfmp*(fmp- f1) + coyfpp)
             / (cstfa2*(f2 - f1));
        y1 =   (yfmp*(f2-fmp) - coyfpp)
             / (cstfa1*(f2 - f1));
      }
      else {
        y2 = yfmp;
        y1 = yfmp;
      }
      if ((fmp-yfmp) < epsi) {
        y2 = f2;
        y1 = f1;
      }
      if (yfmp - std::max(0., (fmp -fs[0])/(1. - fs[0])) < epsi) {
        y2 = std::max(0., (f2 -fs[0]) / (1. - fs[0]));
        y1 = std::max(0., (f1 -fs[0]) / (1. - fs[0]));
      }

      // Determine min and max values of variances of y at F1 and F2.

      y2pmax[0] = (y1-ymin[0])*(ymax[0] - y1);
      y2pmin[0] = 0.;
      y2pmax[1] = (y2-ymin[1])*(ymax[1] - y2);
      y2pmin[1] = 0.;

      // Compute variance at Y max with Y2PMAX 1 and 2

      cs_real_t yfp2max =   cstfa1*(cs_math_pow2(y1) + y2pmax[0])
                          + cstfa2*(cs_math_pow2(y2) + y2pmax[1])
                          - cs_math_pow2(yfmp);
      cs_real_t vymx = yfp2max;

      // Ratio of conditional variances

      cs_real_t cstvar = 0.0;

      if (   ((ymax[1]-y2) > epsi) && ((y2-ymin[1]) > epsi)
          && ((ymax[0]-y1) > epsi) && ((y1-ymin[0]) > epsi)) {
        cstvar =   ((ymax[1]-y2)*(y2-ymin[1]))
                 / ((ymax[0]-y1)*(y1-ymin[0]));
      }

      // Clip VARIANCE Y
      // (we can either clip variance at based on extrema of conditional
      // variances or directly clip conditional variances).

      if (yfp2mp > vymx) {
        if ((yfp2mp-vymx) > mxcyfp) {
          mxcyfp = (yfp2mp-vymx);
          mxyfp2 = vymx;
        }
        yfp2mp = vymx;
        clyf21 += 1;
      }

      // Compute conditional variances

      if (   ((ymax[1]-y2) > epsi) && ((y2-ymin[1]) > epsi)
          && ((ymax[0]-y1) > epsi) && ((y1-ymin[0]) > epsi)) {
        y2p[0] =   (  cs_math_pow2(yfmp) + yfp2mp
                    - cstfa2*cs_math_pow2(y2)-cstfa1*cs_math_pow2(y1))
                 / (cstfa1 + cstfa2*cstvar);
        y2p[1] =   (  cs_math_pow2(yfmp) + yfp2mp
                    - cstfa2*cs_math_pow2(y2) - cstfa1*cs_math_pow2(y1))
                 / (cstfa1/cstvar + cstfa2);
      }
      else if ((ymax[1]-y2 > epsi) && (y2-ymin[1] > epsi)) {
        y2p[0] = 0.;
        y2p[1] =   (  (cs_math_pow2(yfmp)) + yfp2mp
                    - cstfa2*cs_math_pow2(y2) - cstfa1*cs_math_pow2(y1))
                 / cstfa2;
      }
      else if ((ymax[0]-y1 > epsi) && (y1-ymin[0] > epsi)) {
        y2p[1] = 0.;
        y2p[0] =   (  (cs_math_pow2(yfmp)) + yfp2mp
                    - cstfa2*cs_math_pow2(y2) - cstfa1*cs_math_pow2(y1))
                 / cstfa1;
      }
      else {
        y2p[0] = 0.;
        y2p[1] = 0.;
      }

      // Clipping for conditional variances.

      if (y2p[0] > y2pmax[0]) {
        if ((y2p[0] - y2pmax[0]) > mcy2p1) {
          mcy2p1 = (y2p[0]-y2pmax[0]);
          my2p1 = y2pmax[0];
        }
        y2p[0] = y2pmax[0];
        y2p[1] =   (((cs_math_pow2(yfmp)) + yfp2mp
                     -cstfa1*(cs_math_pow2(y1)+y2p[0]))/cstfa2)
                 - cs_math_pow2(y2);
        cly2p1 += 1;
      }
      else if (y2p[0] < y2pmin[0]) {
        if ((y2pmin[0]-y2p[0]) > mcy2p3) {
          mcy2p3 = y2pmin[0] - y2p[0];
          my2p3=y2pmin[0];
        }
        y2p[0] = y2pmin[0];
        y2p[1] =   (((cs_math_pow2(yfmp)) + yfp2mp
                     -cstfa1*(cs_math_pow2(y1)+y2p[0]))/cstfa2)
          - cs_math_pow2(y2);
        cly2p3 += 1;
      }
      if (y2p[1] > y2pmax[1]) {
        if ((y2p[1]-y2pmax[1]) > mcy2p2) {
          mcy2p2 = y2p[1]-y2pmax[1];
          my2p2 = y2pmax[1];
        }
        y2p[1] = y2pmax[1];
        y2p[0] =  (((cs_math_pow2(yfmp)) + yfp2mp
                    -cstfa2*(cs_math_pow2(y2)+y2p[1]))/cstfa1)
                 - cs_math_pow2(y1);
        cly2p2 += 1;
      }
      else if (y2p[1] < y2pmin[1]) {
        if ((y2pmin[1]-y2p[1]) > mcy2p4) {
          mcy2p4 = y2pmin[1] - y2p[1];
          my2p4 = y2pmin[1];
        }
        y2p[1] = y2pmin[1];
        y2p[0] =   (((cs_math_pow2(yfmp))+yfp2mp
                     -cstfa2*(cs_math_pow2(y2)+y2p[1]))/cstfa1)
                 - cs_math_pow2(y1);
        cly2p4 += 1;
      }

      // Compute positions and amplitudes of peaks at Y on F1 and F2.

      _lwcurl(cstfa1, y1, y2p[0], ymin[0], ymax[0],
              y[0], y[1], d[0], d[1]);

      _lwcurl(cstfa2, y2, y2p[1], ymin[1], ymax[1],
              y[2], y[3], d[2], d[3]);

      /* Determine thermochemical quantities of both peaks
         ------------------------------------------------- */

      // Compute enthalpies at 1 and 2.

      cs_real_t sum1  = 0., sum2  = 0., sum3  = 0., sum4 = 0.;
      cs_real_t sum5 = 0., sum6 = 0., sum16 = 0.;

      // Loop on each Dirac.

      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {

        // Compute Enthalpye

        h[dirac_id] = (  (h_max-h_min) * f[dirac_id]
                       + h_min*f_max - h_max*f_min) / (f_max-f_min);

        // Compute mass fraction of gases (F, O and P) for each peak.

        cs_real_t yfuel = y[dirac_id];
        cs_real_t yoxyd = 1. - (coeff3+1.)*f[dirac_id] + coeff3*y[dirac_id];
        cs_real_t yprod = 1. - yfuel - yoxyd;
        cs_real_t yo2   =   coeff1 - (coeff1 + coeff2) * f[dirac_id]
                          + coeff2 * y[dirac_id];

        // Compute molar mass and temperature for each peak.

        coefg[0] = yfuel;
        coefg[1] = yoxyd;
        coefg[2] = yprod;

        // Molar mass for each peak.

        cs_real_t nbmol = 0;
        for (int igg = 0; igg < n_gas_species; igg++)
          nbmol += coefg[igg]/wmolg[igg];
        maml[dirac_id] = 1./nbmol;

        // Temperature for each peak
        teml[dirac_id] = cs_gas_combustion_h_to_t(coefg, h[dirac_id]);

        // Density for each peak
        if (update_rho)
          rhol[dirac_id] =   p0 * maml[dirac_id]
                           / (cs_physical_constants_r*teml[dirac_id]);
        else
          rhol[dirac_id] = ro0;

        // Compute source term for scalar YFM for each peak.
        theta[dirac_id] =   ta / teml[dirac_id]
                          * (1. - teml[dirac_id] / tstar);
        tetmax = std::max(theta[dirac_id], tetmax);
        tetmin = std::min(theta[dirac_id], tetmin);
        w[dirac_id] = vref / lref * (- d[dirac_id]
                                     * yfuel*yo2
                                     * exp(-theta[dirac_id]));

        // Compute source term of scalar YFM for each peak.

        w[dirac_id] =   vref / lref * (- d[dirac_id]
                                       * yfuel*yo2
                                       * exp(-theta[dirac_id]));

        wmax = std::max(w[dirac_id], wmax);
        wmin = std::min(w[dirac_id], wmin);
        o2max = std::max(yo2, o2max);
        o2min = std::min(yo2, o2min);

        // Control sign of W
        w[dirac_id] = std::min(w[dirac_id], 0.);

        // Molar mass of mixture
        sum1 += d[dirac_id]*maml[dirac_id];

        // Mixture temperature
        sum2 += d[dirac_id]*teml[dirac_id];

        // Temperature / Molar mass

        sum3 += d[dirac_id]*teml[dirac_id]/maml[dirac_id];

        // Mass fractions of global species

        sum4 += yfuel*d[dirac_id];
        sum5 += yoxyd*d[dirac_id];
        sum6 += yprod*d[dirac_id];
        sum16 += w[dirac_id];

        // Store properties

        cpro_ampl[dirac_id][c_id] = d[dirac_id];
        cpro_fmel[dirac_id][c_id] = f[dirac_id];
        cpro_fmal[dirac_id][c_id] = y[dirac_id];
        cpro_maml[dirac_id][c_id] = maml[dirac_id];
        cpro_teml[dirac_id][c_id] = teml[dirac_id];
        cpro_rhol[dirac_id][c_id] = rhol[dirac_id];
        cpro_tscl[dirac_id][c_id] = w[dirac_id];

      } // Loop on Diracs

      cpro_mam[c_id]  = sum1;
      cpro_temp[c_id] = sum2;
      cpro_ym1[c_id]  = sum4;
      cpro_ym2[c_id]  = sum5;
      cpro_ym3[c_id]  = sum6;
      cpro_tsc[c_id]  = sum16;

      // Mixture density

      if (update_rho) {
        cs_real_t temsmm = sum3;
        crom[c_id] =    srrom * crom[c_id]
                     + (1.-srrom) * (p0/(cs_physical_constants_r*temsmm));
      }

    } // End if PDF case

  } // End of loop on cells

  /* Logging
     ------- */

  if (cs_log_default_is_active() == false)
    return;

  cs_parall_sum_scalars(clfp2m, clyf21, clcyf1, clcyf2, cly2p1, cly2p3,
                        cly2p2, cly2p4, icpt1, icpt2);
  cs_parall_max_scalars(mxcfp, maxfp2, mxcyfp, mxyfp2, mxccyf, mxcoyf,
                        mnccyf, mcy2p1, my2p1, mcy2p3, mcy2p2, my2p2,
                        mcy2p4, o2max, tetmax, wmax);
  cs_parall_min_scalars(mncoyf, my2p3, my2p4, o2min, tetmin, wmin);

  cs_log_printf(CS_LOG_DEFAULT,
                _(" clips to high on variance at f = %lu\n"
                  " maximum extent (value reached - max value) = %g\n"
                  " max value (for maximum extent) = %g\n"
                  "\n"
                  " clips to high on variance at y = %lu\n"
                  " maximum extent (value reached - max value) = %g\n"
                  " max value (for maximum extent) = %g\n"
                  "\n"
                  " clips to high on covariance = %lu\n"
                  " maximum extent (value reached - max value) F = %g\n"
                  " max value (for maximum extent) = %g\n"
                  "\n"
                  " clips to low on covariance = %lu\n"
                  " maximum extent (value reached - min value) F = %g\n"
                  " min value (for maximum extent) = %g\n"
                  "\n"),
                (unsigned long)clfp2m, mxcfp, maxfp2,
                (unsigned long)clyf21, mxcyfp, mxyfp2,
                (unsigned long)clcyf1, mxccyf, mxcoyf,
                (unsigned long)clcyf2, mnccyf, mncoyf);

  cs_log_printf(CS_LOG_DEFAULT,
                _(" clips to high on conditional variance 1 = %lu\n"
                  " maximum extent (value reached - max value) = %g\n"
                  " max value (for maximum extent) = %g\n"
                  "\n"
                  " clips to low on conditional variance 1 = %lu\n"
                  " maximum extent (value reached - min value) = %g\n"
                  " min value (for maximum extent) = %g\n"
                  "\n"
                  " clips to high on conditional variance 2 = %lu\n"
                  " maximum extent (value reached - max value) = %g\n"
                  " max value (for maximum extent) = %g\n"
                  "\n"
                  " clips to low on conditional variance 2 = %lu\n"
                  " maximum extent (value reached - min value) = %g\n"
                  " min value (for maximum extent) = %g\n"
                  "\n"),
                (unsigned long)cly2p1, mcy2p1, my2p1,
                (unsigned long)cly2p3, mcy2p3, my2p3,
                (unsigned long)cly2p2, mcy2p2, my2p2,
                (unsigned long)cly2p4, mcy2p4, my2p4);

  cs_log_printf(CS_LOG_DEFAULT,
                _(" points without PDF = %lu\n"
                  " points with    PDF = %lu\n"
                  " min max 02    = %g %g\n"
                  " min max theta = %g %g\n"
                  " min max W     = %g %g\n"
                  "\n"),
                (unsigned long)icpt1, (unsigned long)icpt2,
                o2min, o2max, tetmin, tetmax, wmin, wmax);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for Libby-Williams gas combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_fields_init(void)
{
  // Only when not a restart
  if (cs_restart_present())
    return;

  // Local variables

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int igg = 0; igg < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; igg++)
    coefg[igg] = 0;

  // Initializations with air at tinitk
  // ----------------------------------

  // Mixture enthalpy
  if (cm->type % 2 == 1) {
    // mixture temperature: air at tinitk
    cs_real_t tinitk = cs_glob_fluid_properties->t0;

    // Air enthalpy at tinitk
    coefg[0] = 0.;
    coefg[1] = 1.;
    coefg[2] = 0.;
    cs_real_t hair = cs_gas_combustion_t_to_h(coefg, tinitk);

    // Mixture enthalpy
    cs_real_t *cvar_scalt = CS_F_(h)->val;
    cs_array_real_set_scalar(n_cells_ext, hair, cvar_scalt);
  }

  // Mass fraction
  cs_array_real_set_scalar(n_cells_ext, cm->lw.fmax, cm->yfm->val);

  // Mixture fraction
  cs_array_real_set_scalar(n_cells_ext, cm->lw.fmax, cm->fm->val);

  // No need to set yfp2m, fp2m, or coyfp to 0,
  // as this is the default for all fields.
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute physical properties for Libby-Williams combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_physical_prop(void)
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const int sub_type = cm->type % 100;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_t *cvar_fm = cm->fm->val;
  cs_real_t *cvar_fp2m = cm->fp2m->val;
  cs_real_t *cvar_yfm = cm->yfm->val;
  cs_real_t *cvar_yfp2m = cm->yfp2m->val;

  cs_real_t *cvar_coyfp = nullptr;
  if (cm->coyfp != nullptr)
    cvar_coyfp = cm->coyfp->val;

  /*
   * Determine thermochemical quantities
   * ----------------------------------- */

  if (sub_type == 0 || sub_type == 1)
    _pdflwc(n_cells, cvar_fm, cvar_fp2m, cvar_yfm, cvar_yfp2m);

  else if (sub_type == 2 || sub_type == 3)
    _pdfpp3(n_cells, cvar_fm, cvar_fp2m, cvar_yfm, cvar_yfp2m, cvar_coyfp);

  else if (sub_type == 4 || sub_type == 5)
    _pdfpp4(n_cells, cvar_fm, cvar_fp2m, cvar_yfm, cvar_yfp2m, cvar_coyfp);

  /*
   * Compute rho and mass fractions of global species at boundaries
   * -------------------------------------------------------------- */

  cs_combustion_boundary_conditions_density_ebu_lw();

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  cs_host_context ctx;

  for (int igg = 0; igg < cm->n_gas_species; igg++) {
    cs_real_t *bsval = cm->bym[igg]->val;
    cs_real_t *cpro_ymgg = cm->ym[igg]->val;

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST (cs_lnum_t f_id) {
      cs_lnum_t c_id = b_face_cells[f_id];
      bsval[f_id] = cpro_ymgg[c_id];
    });
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute source terms for premixed flame Libby-Williams
 *        combustion model.
 *
 * Define the source terms for a given scalar over one time step.
 *
 * The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 * \f$ rovsdt \f$ et \f$ smbrs \f$ may already contain source term
 * so must not be overwritten, but incremented.
 *
 * For stability, only positive terms should be add in \f$ rovsdt \f$.
 * There is no constraint for \f$ smbrs \f$.
 * For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
 *           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
 *           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
 *
 * Here we set \f$ rovsdt \f$ and \f$ smbrs \f$ containing \f$ \rho \Omega \f$
 *   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
 *     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.C.s^{-1} \f$,
 *     for enthalpy: \f$ J.s^{-1} \f$)
 *   - \f$ rovsdt \f$ in \f$ kg.s^{-1} \f$
 *
 * \param[in]      f_sc          pointer to scalar field
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_source_terms(cs_field_t  *f_sc,
                              cs_real_t    smbrs[],
                              cs_real_t    rovsdt[])
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_sc);
  if (eqp->verbosity >= 1)
    bft_printf(_("Source terms for variable %s\n\n"), f_sc->name);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  // Get variables and coefficients

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t *cvara_k = nullptr;
  const cs_real_t *cvara_ep = nullptr;
  const cs_real_t *cvara_omg = nullptr;
  const cs_real_6_t *cvara_rij = nullptr;

  const cs_real_t *cvara_scal = f_sc->val_pre;
  const cs_real_t *cvara_yfm = cm->yfm->val_pre;
  const cs_real_t *cvara_fm = cm->fm->val_pre;

  cs_real_t *cpro_fmel[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_fmal[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_tscl[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_rhol[CS_COMBUSTION_GAS_MAX_DIRAC];

  const int n_dirac = cm->lw.n_dirac;
  for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
    cpro_fmel[dirac_id] = cm->lw.fmel[dirac_id]->val;
    cpro_fmal[dirac_id] = cm->lw.fmal[dirac_id]->val;
    cpro_rhol[dirac_id] = cm->lw.rhol[dirac_id]->val;
    cpro_tscl[dirac_id] = cm->lw.tscl[dirac_id]->val;
  }

  const cs_real_t epsi = 1.0e-10;

  cs_host_context ctx;

  /* Source term for mean fuel mass fraction
     --------------------------------------- */

  if (f_sc == cm->yfm) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t sum = 0.;
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
        sum +=   cpro_rhol[dirac_id][c_id]
               * cpro_tscl[dirac_id][c_id] * volume[c_id];
      }

      // implicit term
      if (cvara_scal[c_id] > epsi)
        rovsdt[c_id] += std::max(-sum/cvara_scal[c_id], 0.0);

      // explicit term
      smbrs[c_id] += sum;
    });

  }

  /* Source term for variance of mean fuel mass fraction
     --------------------------------------------------- */

  if (f_sc == cm->yfp2m) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t sum = 0.;
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
        sum  +=  (cpro_tscl[dirac_id][c_id]*volume[c_id]
                *(cpro_fmal[dirac_id][c_id] - cvara_yfm[c_id])
                *cpro_rhol[dirac_id][c_id]);
      }
      smbrs[c_id] += sum;
    });

  }

  /* Covariance source term
     ---------------------- */

  if (f_sc == cm->coyfp) {

    // Allocate work arrays

    cs_real_3_t *gradf, *grady;
    CS_MALLOC(gradf, n_cells_ext, cs_real_3_t);
    CS_MALLOC(grady, n_cells_ext, cs_real_3_t);

    // Gradient of F
    cs_field_gradient_scalar(cm->fm, true, 1, gradf);

    // Gradient of Yfuel
    cs_field_gradient_scalar(cm->yfm, true, 1, grady);

    // Access or reconstruct k and Epsilon based on turbulence model.

    cs_real_t  *w1 = nullptr;

    if (   cs_glob_turb_model->itytur == 2
        || cs_glob_turb_model->itytur == 5) {
      cvara_k = CS_F_(k)->val_pre;
      cvara_ep = CS_F_(eps)->val_pre;
    }

    else if (cs_glob_turb_model->itytur == 3) {
      cvara_rij = (const cs_real_6_t *)(CS_F_(k)->val_pre);
      cvara_ep = CS_F_(eps)->val_pre;

      CS_MALLOC(w1, n_cells_ext, cs_real_t);
      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        w1[c_id] = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
      });
      cvara_k = w1;
    }

    else if (cs_glob_turb_model->model == CS_TURB_K_OMEGA) {
      cvara_k = CS_F_(k)->val_pre;
      cvara_omg = CS_F_(omg)->val_pre;

      const cs_real_t cmu = cs_turb_cmu;
      CS_MALLOC(w1, n_cells_ext, cs_real_t);
      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        w1[c_id] = cmu * cvara_k[c_id]* cvara_omg[c_id];
      });
      cvara_ep = w1;
    }

    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const int krvarfl = cs_field_key_id("variance_dissipation");
    cs_real_t turb_schmidt = cs_field_get_key_double(f_sc, ksigmas);
    cs_real_t rvarfl = cs_field_get_key_double(f_sc, krvarfl);

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

      // To confirm:
      // The dissipation term should be implicit.
      // In the dissipation term, a constant Cf is missing.
      //   Can it be conidered == 1 ?
      // Check the sign of the production term.

      // Implicit term

      cs_real_t w11 =   cvara_ep[c_id] / (cvara_k[c_id] * rvarfl)
                      * volume[c_id] * crom[c_id];
      rovsdt[c_id] = rovsdt[c_id] + std::max(w11, 0.);

      // Gradient term

      cs_real_t tsgrad =  (  2.0
                           * visct[c_id]/(turb_schmidt)
                           * cs_math_3_dot_product(gradf[c_id], grady[c_id]))
                         * volume[c_id];

      // Dissipation term

      cs_real_t tsdiss = -w11 * cvara_scal[c_id];

      // Chemical term

      cs_real_t tschim = 0.;
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
        tschim += ( cpro_tscl[dirac_id][c_id]
                   *(cpro_fmel[dirac_id][c_id]-cvara_fm[c_id])
                   *volume[c_id]) * cpro_rhol[dirac_id][c_id];
      }

      // Sum of sources

      smbrs[c_id] += tschim + tsgrad + tsdiss;
    });

    CS_FREE(w1);
    CS_FREE(grady);
    CS_FREE(gradf);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
