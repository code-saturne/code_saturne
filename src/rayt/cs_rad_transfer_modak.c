/*============================================================================
 * Radiation solver MODAK library.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_parameters.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_modak.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_rad_transfer_modak.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for fortran API
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute Chebychev polynomial of order norpol
 *
 * \param[in]   norpol    chebychev polynomial order
 * \param[in]   argpol    chebychev polynomial argument
 *
 * \returns     valpol    chebychev polynomial value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_chebyc(int        norpol,
        cs_real_t  argpol)
{
  cs_real_t valpol;
  if (norpol <= 0)
    valpol = 1.0;
  else if (norpol <= 1)
    valpol = argpol;
  else {
    cs_real_t f = argpol + argpol;
    cs_real_t vm1 = argpol;
    cs_real_t vm2 = 1.0;
    for (int ict = 1; ict < norpol; ict++) {
      valpol = f * vm1 - vm2;
      vm2 = vm1;
      vm1 = valpol;
    }
  }
  return valpol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute asymptotic expression for function pentagamma
 *
 * \param[in] zz
 *
 * \return    Pentagamma function value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_asympt(cs_real_t  zz)
{
  cs_real_t d1s3 = 1.0 / 3.0;
  cs_real_t zi1 = 1.0 / zz;
  cs_real_t zi2 = zi1 * zi1;
  cs_real_t zi3 = zi1 * zi2;
  return   zi3
         * (  (2.0 + 3.0 * zi1)
            + zi2 * (2.0 + zi2 * (- 1.0
                                  + zi2 * (  1.0 + d1s3
                                           + zi2 * ( -3.0 + 10.0 * zi2))
                                 )
                    )
           );
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute value of pentagamma function of argument x.
 *
 * We use asymptotic and recurrence formulas of Abramowitz and Stegun.
 *
 * \param[in]  argfpe   pentagamma function argument
 *
 * \return  pentagamma function value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_pentag(cs_real_t  argfpe)
{
  cs_real_t  zz, zzv, zs;
  if (argfpe >= 4.0) {
    zs = 0.0;
    zz = argfpe;
  }
  else if (argfpe >= 3.0) {
    zs = 6.0 / pow(argfpe, 4.0);
    zz = argfpe + 1.0;
  }
  else if (argfpe >= 2.0) {
    zs = (  1.0 / pow (argfpe + 1.0, 4.0)
          + 1.0 / pow (argfpe      , 4.0)) * 6.0;
    zz = argfpe + 2.0;
  }
  else {
    zs  = (  1.0 / pow(argfpe + 2.0, 4.0)
           + 1.0 / pow(argfpe + 1.0, 4.0)
           + 1.0 / pow(argfpe      , 4.0)) * 6.0;
    zz  = argfpe + 3.0;
  }

  zzv = _asympt(zz);

  return zzv + zs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute transmissivity (taus) of path at a given temperature.
 *
 *\param[in]  zkled
 *\param[in]  pathl    penetration of radiation in mixture
 *\param[in]  tblack   source or gas temperature
 *
 *\return taus     transmissivity of path
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_tasoot(cs_real_t  zkled,
        cs_real_t  pathl,
        cs_real_t  tblack)
{
  if (zkled <= 0.0)
    return 1.0;
  else {
    cs_real_t arg = 1.0 + zkled * pathl * tblack * 6.5333e-05;
    return _pentag(arg) * 0.1539897336;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute correction for mi of co2 et h2o when wavelengths
 *        are beyond 2.7 and 15 microns.
 *
 * \param[in] val
 * \param[in] pl
 * \param[in] te
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_fdleck(cs_real_t  val,
        cs_real_t  pl,
        cs_real_t  te)
{
  cs_real_t term, term2, term3, tt, tt2, aa, bb, cc;

  if (pl >= 0.1) {

    term   = val / (10.7 + 101.0 * val) - pow(val, 10.4) / 111.7;
    term2  = pow(log10 (101.325 * pl), 2.76);
    tt     = te / 1000.0;
    tt2    = tt * tt;
    aa     =  -1.0204082;
    bb     = 2.2448979;
    cc     =  -0.23469386;
    term3  = aa * tt2 + bb * tt + cc;

    /* term3 represents the temperature adjustment */

    return term * term2 * term3;
  }
  else
    return 0.0;
}


/*----------------------------------------------------------------------------*/
/*!
 * \param[in]  pp
 * \param[in]  pl
 * \param[in]  te
 * \param[in]  index
 * \param[out] val
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_scrtch(cs_real_t  pp,
        cs_real_t  pl,
        cs_real_t  te,
        int        index)
{
  cs_real_t cc[3][4][4], cw[3][4][4];

  /* CC represents an array of 48 elements for CO2 */
  if (index != 2) {
    cc[0][0][0] = -2.754568;
    cc[0][0][1] = -0.2997857;
    cc[0][0][2] = -0.1232494;
    cc[0][0][3] =  0.01279287;
    cc[0][1][0] =  1.503051;
    cc[0][1][1] =  0.3156449;
    cc[0][1][2] =  0.01058126;
    cc[0][1][3] = -0.03729625;
    cc[0][2][0] = -0.247411;
    cc[0][2][1] = -0.03323846;
    cc[0][2][2] = -0.01819471;
    cc[0][2][3] =  0.02289789;
    cc[0][3][0] =  0.04994029;
    cc[0][3][1] = -0.001986786;
    cc[0][3][2] =  0.003007898;
    cc[0][3][3] = -0.001175598;
    cc[1][0][0] =  0.005737722;
    cc[1][0][1] = -0.009328458;
    cc[1][0][2] =  0.002906286;
    cc[1][0][3] =  0.000422752;
    cc[1][1][0] = -0.003151784;
    cc[1][1][1] =  0.005632821;
    cc[1][1][2] = -0.003260295;
    cc[1][1][3] =  0.0007065884;
    cc[1][2][0] =  0.0001668751;
    cc[1][2][1] = -0.0007326533;
    cc[1][2][2] =  0.0003639855;
    cc[1][2][3] =  0.0003228318;
    cc[1][3][0] =  0.0007386638;
    cc[1][3][1] = -0.0007277073;
    cc[1][3][2] =  0.0005925968;
    cc[1][3][3] = -0.0002021413;
    cc[2][0][0] =  0.003385611;
    cc[2][0][1] = -0.005439185;
    cc[2][0][2] =  0.00176456;
    cc[2][0][3] =  0.0003036031;
    cc[2][1][0] = -0.0018627;
    cc[2][1][1] =  0.003236275;
    cc[2][1][2] = -0.00195225;
    cc[2][1][3] =  0.0003474022;
    cc[2][2][0] =  0.0001204807;
    cc[2][2][1] = -0.0004479927;
    cc[2][2][2] =  0.0002497521;
    cc[2][2][3] =  0.0001812996;
    cc[2][3][0] =  0.0004218169;
    cc[2][3][1] = -0.0004046608;
    cc[2][3][2] =  0.0003256861;
    cc[2][3][3] = -9.514981e-05;
  }

  /* CW represents an array of 48 elements for H2O */
  else {
    cw[0][0][0] = -2.594279;
    cw[0][0][1] = -0.7118472;
    cw[0][0][2] = -0.0009956839;
    cw[0][0][3] =  0.0122656;
    cw[0][1][0] =  2.510331;
    cw[0][1][1] =  0.6481808;
    cw[0][1][2] = -0.03330587;
    cw[0][1][3] = -0.005524345;
    cw[0][2][0] = -0.4191636;
    cw[0][2][1] = -0.137518;
    cw[0][2][2] =  0.0387793;
    cw[0][2][3] =  0.0008862328;
    cw[0][3][0] = -0.0322912;
    cw[0][3][1] = -0.01820241;
    cw[0][3][2] = -0.02223133;
    cw[0][3][3] = -0.0005940781;
    cw[1][0][0] =  0.1126869;
    cw[1][0][1] = -0.08133829;
    cw[1][0][2] =  0.0151494;
    cw[1][0][3] =  0.00139398;
    cw[1][1][0] = -0.009298805;
    cw[1][1][1] =  0.0455066;
    cw[1][1][2] = -0.02082008;
    cw[1][1][3] =  0.002013361;
    cw[1][2][0] = -0.04375032;
    cw[1][2][1] =  0.01924597;
    cw[1][2][2] =  0.008859877;
    cw[1][2][3] = -0.004618414;
    cw[1][3][0] =  0.007077876;
    cw[1][3][1] = -0.02096188;
    cw[1][3][2] =  0.001458262;
    cw[1][3][3] =  0.003851421;
    cw[2][0][0] =  0.05341517;
    cw[2][0][1] = -0.03407693;
    cw[2][0][2] =  0.004354611;
    cw[2][0][3] =  0.001492038;
    cw[2][1][0] = -0.004708178;
    cw[2][1][1] =  0.02086896;
    cw[2][1][2] = -0.009477533;
    cw[2][1][3] =  0.0006153272;
    cw[2][2][0] = -0.02104622;
    cw[2][2][1] =  0.007515796;
    cw[2][2][2] =  0.005965509;
    cw[2][2][3] = -0.002756144;
    cw[2][3][0] =  0.004318975;
    cw[2][3][1] = -0.01005744;
    cw[2][3][2] =  0.0004091084;
    cw[2][3][3] =  0.002550435;
  }

  cs_real_t xx = log(pp) / 3.45 + 1.0;
  cs_real_t yy = (log(pl) + 2.555) / 4.345;
  cs_real_t zz = (te - 1150.0) / 850.0;

  cs_real_t value = 0.0;

  for (int ii = 0; ii < 3; ii++) {

    cs_real_t tix;
    tix = _chebyc(ii, xx);

    cs_real_t v6 = 0.0;
    for (int jj = 0; jj < 4; jj++) {

      cs_real_t tjy;
      tjy = _chebyc(jj, yy);
      cs_real_t v7 = 0.0;

      for (int kk = 0; kk < 4; kk++) {

        cs_real_t tkz = _chebyc(kk, zz);

        if (index == 1)
          v7 += tkz * cc[ii][jj][kk];
        if (index == 2)
          v7 += tkz * cw[ii][jj][kk];


      }
      v6 += v7 * tjy;
    }

    value += v6 * tix;
  }

  return exp(value);
}

/*----------------------------------------------------------------------------*/
/*!
 *\brief Compute emmissivity for given path of a mix of co2 ans
 *       h2o at temperature te
 *
 * \param[in]  pathl value of path
 * \param[in]  pc    co2 partial pressure
 * \param[in]  pw    h2o partial pressure
 * \param[in]  te    temperature
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_emigas(cs_real_t  pathl,
        cs_real_t  pc,
        cs_real_t  pw,
        cs_real_t te)
{
  cs_real_t tmin   = 298.0;
  cs_real_t tmax   = 3000.0;
  cs_real_t val_emigas  = 0.0;
  if (te < tmin || te > tmax)
    return 0.0;

  cs_real_t ec = 0.0;
  if (pc >= 0.0011 && pc <= 1.0) {

    cs_real_t pcl = pc * pathl;

    if ( pcl >= 0.0011 && pcl <= 5.98)
      ec = _scrtch(pc, pcl, te, 1);

  }

  if (pw >= 0.0011 && pw <= 1.0) {

    cs_real_t pwl = pw * pathl;

    if (pwl >= 0.0011 && pwl <= 5.98) {

      cs_real_t ew = _scrtch (pw, pwl, te, 2);
      val_emigas = ec + ew;

      if (ec <= 0.0)
        return val_emigas;

      cs_real_t pcpw = pc + pw;
      cs_real_t xi   = pw / pcpw;
      if (xi < 0.01)
        return val_emigas;

      cs_real_t pcwl = pcpw * pathl;
      if (pcwl < 0.1)
        return val_emigas;

      cs_real_t dels = _fdleck(xi,pcwl,te);
      val_emigas -= dels;
      return val_emigas;
    }

  }
  val_emigas  = ec;

  return val_emigas;
}

/*----------------------------------------------------------------------------*/
/*!
 *   on calcule les absorptivites (par rapport
 *   a une source de corps noir) d'un melange gazeux isotherme,
 *   homogene de suie, co2 et h2o a la pression totale d'1 atm.
 *   si la temperature du corps noir est egale a la temperature du
 *   melange, l'absorptivite est egale a l'emissivite.
 *   les emissivites ainsi calculees sont en bon accord avec les
 *   calculs spectraux et les mesures experimentales
 *   ts et te doivent etre compris entre 300 et 2000 kelvin
 *   la longueur d'onde 0.94 micron.
 *   sootk est liee a la fraction volumique de suie fv suivant
 *   la formule :
 *                  SOOTK=7FV/0.94E-6
 *
 * \param[in]   ts       temperature du corps noir (k)
 * \param[in]   te       temperature du melange (k)
 * \param[in]   path     penetration du rayonnement dans le melange (m)
 * \param[in]   sootk    coefficient d'absorption des suies
 * \param[in]   pco2     pression partielle de co2 dans un melange
 *                                       de presssion totale 1 atm.
 * \param[in]   ph2o     pression partielle de h2o dans un melange
 *                                       de presssion totale 1 atm.
 * \param[in]   alpha    absorptivite
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
_absorb(cs_real_t  ts,
        cs_real_t  te,
        cs_real_t  path,
        cs_real_t  sootk,
        cs_real_t  pco2,
        cs_real_t  ph2o)
{
  cs_real_t tmax   = 3000.0;
  cs_real_t tmin   = 298.0;

  if (   ts > tmin && ts < tmax
      && te > tmin && te < tmax) {

    /* --- Pression totale : PTOTAL   */
    cs_real_t ptotal = pco2 + ph2o;

    /* --- Rapport temeperature melange et temperature source : RATIO    */
    if (ptotal <= 1.0) {

      cs_real_t ratio = te / ts;

      /* --- Longueur de penetration du rayonment effectif : PATHL    */
      cs_real_t pathl = path / ratio;
      cs_real_t pcl   = pco2 * pathl;
      cs_real_t pwl   = ph2o * pathl;

      if ( pcl < 5.98 && pwl < 5.98) {

        /* --- Calcul de l'absortivite des suies : AS    */
        cs_real_t as  = 0.0;
        if (sootk > 0.0)
          as = 1.0 - _tasoot(sootk, path, ts);

        /* --- Calcul de l'absorptivite du gaz : AG */
        /*                 = emissivite du gaz */
        cs_real_t ag  = 0.0;
        if (   (pco2 >= 0.0011 || ph2o >= 0.0011)
            && (pcl >= 0.0011 || pwl >= 0.0011)) {
          ag = _emigas(pathl, pco2 ,ph2o, ts);

          /* --- Calcul de la fraction de vapeur d'eau : ZETA   */
          cs_real_t zeta = ph2o / ptotal;

          cs_real_t power = 0.65 - 0.2 * zeta;
          ag = ag * (pow (ratio, power));
        }

        cs_real_t alpha = as + ag - as * ag;

        if (alpha > 1e-08)
          return alpha;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("\n"
                    "Modak model error:\n"
                    "  the product path*Ts/T*pCO2 or path*Ts/T*pH2O\n"
                    "  is greater than 5.98 atm.meters."));
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("\n"
                  "Modak model error:\n"
                  "  the sum of partial pressures of CO2 and H2O gases\n"
                  "  is greater than 1 atmosphere."));
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("\n"
                "Modak model error:\n"
                "  the mixture temperature Te or blackbody temperature Ts\n"
                "  is out of domain validity bounds."));

  return 1e-08;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute radiative properties of a gas based on temperature,
 *        composition of products in co2, h2o and soot
 *        using linear regressions established by Modak.
 *
 * \param[out] ck     coefficient d'absorption du milieu (nul si transparent)
 * \param[in]  pco2   pression partielle de co2
 * \param[in]  pco2   pression partielle de h2o
 * \param[in]  fv     fraction volumique de suies
 * \param[in]  temp   temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_modak(cs_real_t        ck[],
                      const cs_real_t  pco2[],
                      const cs_real_t  ph2o[],
                      const cs_real_t  fv[],
                      const cs_real_t  temp[])
{
  /* Mean penetration length of radiation */
  cs_real_t path = 15.0;

  cs_real_t tmax = 2000.0;
  cs_real_t tmin = 300.0;

  /* Caution: temperatures used by Modak are in Kelvin */

  for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {

    cs_real_t te; /* gas mix temperature */
    cs_real_t ts; /* black body temperature */

    /* Limitation to Tmax = 2000 K and Tmin = 300 K */
    if (temp[cell_id] > tmax) {
      ts = tmax;
      te = tmax;
    }
    else if (temp[cell_id] < tmin) {
      ts = tmin;
      te = tmin;
    } else {
      ts = temp[cell_id];
      te = temp[cell_id];
    }

    /* soot volume fraction */
    cs_real_t sootk = 7.0 * fv[cell_id] / 9.5e-07;

    /* fluid absorptivity */
    cs_real_t alpha = _absorb(ts, te, path, sootk, pco2[cell_id], ph2o[cell_id]);

    /* check */
    if ((1.0 - alpha) <= cs_math_epzero)
      bft_error(__FILE__, __LINE__, 0,
                _("Error in %s: absorptivity computation\n"
                  "  cell_id = %10d\n"
                  "  alpha = %15.7e\n"
                  "  pco2  = %15.7e\n"
                  "  ph2o  = %15.7e\n"
                  "  sootk = %15.7e\n"
                  "  te    = %15.7e\n"
                  "  path  = %15.7e\n"
                  "  fv    = %15.7E\n"),
                __func__,
                cell_id,
                alpha,
                pco2[cell_id],
                ph2o[cell_id],
                sootk,
                te,
                path,
                fv[cell_id]);

    /* Compute absorption coefficient */
    ck[cell_id] = -log(1.0 - alpha) / path;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
