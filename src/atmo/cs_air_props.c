/*============================================================================
 * Specific laws for air properties (temperature, enthalpy)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_math.h"
#include "cs_physical_constants.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_air_props.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main physical constants structure and associated pointer */

static cs_air_fluid_props_t _props = {
  .humidity0 = 0.,
  .cp_a = 0.,
  .cp_v = 0.,
  .cp_l = 0.,
  .hv0 = 0.,
  .rho_l = 0.,
  .lambda_h = 0.,
  .lambda_l = 0.,
  .droplet_diam = 0.
};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Global variables
 *============================================================================*/

cs_air_fluid_props_t *cs_glob_air_props = &_props;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_air_glob_properties_get_pointer(double  **humidity0,
                                   double  **cp_a,
                                   double  **cp_v,
                                   double  **cp_l,
                                   double  **hv0,
                                   double  **rho_l,
                                   double  **lambda_h,
                                   double  **lambda_l,
                                   double  **droplet_diam);

/*============================================================================
 *  Public functions for Fortran API
 *============================================================================*/

void
cs_air_glob_properties_get_pointer(double  **humidity0,
                                   double  **cp_a,
                                   double  **cp_v,
                                   double  **cp_l,
                                   double  **hv0,
                                   double  **rho_l,
                                   double  **lambda_h,
                                   double  **lambda_l,
                                   double  **droplet_diam)
{
  *humidity0    = &(_props.humidity0);
  *cp_a         = &(_props.cp_a);
  *cp_v         = &(_props.cp_v);
  *cp_l         = &(_props.cp_l);
  *hv0          = &(_props.hv0);
  *rho_l        = &(_props.rho_l);
  *lambda_h     = &(_props.lambda_h);
  *lambda_l     = &(_props.lambda_l);
  *droplet_diam = &(_props.droplet_diam );
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation water mass enthalpy
 *
 * \return Liquid water mass enthalpy
 *
 * \param[in]     t_l          water temperature in Celsius degree
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_h_l(cs_real_t  t_l)
{
  cs_real_t h_l = (cs_glob_air_props->cp_l) * t_l;

  return  h_l;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation water vapor mass enthalpy
 *
 * \return water vapour mass enthalpy
 *
 * \param[in]     t_vap          water vapour temperature in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_hvap(cs_real_t  t_vap)
{
  cs_air_fluid_props_t  *ct_prop = cs_glob_air_props;

  cs_real_t h_vap;

  h_vap = ct_prop->hv0 + ct_prop->cp_v * t_vap;

  return  h_vap;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the derivate of the absolute humidity at saturation
 *
 * \return derivative of the humidity of saturated air
 *
 * \param[in]     th            temperature in Celsius degree
 * \param[in]     p0            reference pressure
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_dxsath(cs_real_t  th,
              cs_real_t  p0)
{
  cs_real_t   a1,b1,c1,ps,pv,grpim;
  cs_real_t   dxsath = 0.;

  /* T between -20 and 0 degrees C */
  /* Warning if T less than -20 degrees C */
  if (th >= -20. && th <= 0.) {

    a1 = 6.4147;
    b1 = 22.376;
    c1 = 271.68;
    ps = a1 + (b1 * th)/(c1 + th);
    pv = exp(ps);
    grpim =  b1 *   c1 / pow(c1 + th,2.);

    dxsath = 0.622 * pv * p0 * grpim /pow(p0 - pv,2.);

  }

  /* T between 0 and 40 degrees C */

  else if (th >=0. && th <= 40.) {

    a1 = 6.4147;
    b1 = 17.438;
    c1 = 239.78;
    ps = a1 + (b1 * th)/(c1 + th);
    pv = exp(ps);
    grpim =  b1 *   c1 / pow(c1 + th,2.);

    dxsath =0.622 * pv * p0 * grpim /pow(p0 - pv,2.);

  }

  /* T between 40 and 80 degrees C */

  else if (th >= 40. && th <= 80.) {

    const cs_real_t  T0 = 273.16;
    const cs_real_t  Ax = 8.2969;
    const cs_real_t  Ay = 4.76955;
    const cs_real_t  A0 = 0.78614;
    const cs_real_t  A2 = 5.028;
    const cs_real_t  A3 = 0.000150475;
    const cs_real_t  A4 = 0.00042873;
    cs_real_t  tt, px, px10, py, py10, pspr, pvpr,
               g1, g1pr, g2, g2pr, g3, g3pr, g4, g4pr;

    a1 = 10.7954;

    tt = th / T0 ;
    px = Ax * tt;
    px10 = pow( 10., px );
    py = Ay *tt / (1. + tt );
    py10 = pow( 10., py );
    g1 = a1 * tt / (1. + tt);
    g1pr = a1 / (T0 * pow(1. + tt, 2.) );
    g2 = - A2 * log10(1. + tt );
    g2pr = - A2 /(T0 * log(10.)*(1. + tt));
    g3 = A3 *(1. - 1. /px10);
    g3pr = A3 * Ax * log(10.) / (T0 * px10);
    g4 = A4 *(py10 - 1.);
    g4pr = A4 * Ay * log (10.) * py10 / (T0* pow(1. + tt, 2.));
    ps = A0 + g1 + g2 + g3 + g4;
    pspr = g1pr + g2pr + g3pr + g4pr;
    pv = pow(10., ps) *100.;
    pvpr  = log(10.) * pspr * pv;
    dxsath= p0 * pvpr * 0.622 / pow (p0 - pv, 2);

  }
  else if (th > 80.) {

    dxsath =  0.001;

  }

  return dxsath ;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the Cp of humid air
 *
 * \return specific heat of humid air
 *
 * \param[in]     x             absolute humidity of humid air
 * \param[in]     x_s           absolute humidity of saturated humid air
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_cp_humidair(cs_real_t  x,
                   cs_real_t  x_s)
{
  cs_real_t  cp_h;
  cs_air_fluid_props_t  *ct_prop = cs_glob_air_props;

  if (x <= x_s) {
    cp_h = ct_prop->cp_a + x*ct_prop->cp_v;
  } else {
    cp_h = ct_prop->cp_a + x_s*ct_prop->cp_v + (x-x_s)*ct_prop->cp_l;
  } //FIXME + dX_s * hv ?

  cp_h = cp_h / (1.0+x);

  return cp_h;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the specific enthalpy of humid air
 *
 * \return specific enthalpy of humid air
 *
 * \param[in]     cp_h          Cp of humid air
 * \param[in]     x             absolute humidity of humid air
 * \param[in]     x_s           absolute humidity of saturated humid air
 * \param[in]     t_h           temperature of humid air in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_h_humidair(cs_real_t  cp_h,
                  cs_real_t  x,
                  cs_real_t  x_s,
                  cs_real_t  t_h)
{
  cs_air_fluid_props_t  *ct_prop = cs_glob_air_props;
  cs_real_t h_h;

  cs_real_t tk_h = t_h + cs_physical_constants_celsius_to_kelvin;

  if (x <= x_s) {
    h_h  = cp_h*tk_h + x*ct_prop->hv0/(1.0+x);
  } else {
    h_h  = cp_h*tk_h + x_s*ct_prop->hv0/(1.0+x);
  }

  return  h_h;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the temperature of humid air
 *
 * \return temperature of humid air (in Celsius)
 *
 * \param[in]     cp_h          Cp of humid air
 * \param[in]     x             absolute humidity of humid air
 * \param[in]     x_s           absolute humidity of saturated humid air
 * \param[in]     h_h           humid air enthalpy
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_t_humidair(cs_real_t  cp_h,
                  cs_real_t  x,
                  cs_real_t  x_s,
                  cs_real_t  h_h)
{
  cs_air_fluid_props_t  *ct_prop = cs_glob_air_props;
  cs_real_t t_h;

  if (x <= x_s) {
    t_h = (h_h - x     * ct_prop->hv0/(1.0+x))/cp_h;
  } else {
    t_h = (h_h - x_s * ct_prop->hv0/(1.0+x))/cp_h;
  }

  return  t_h;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the temperature of liquid water
 *
 * \return liquid water temperature (in Celsius)
 *
 * \param[in]     h_l           specific enthalpy of liquid water
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_liq_h_to_t(cs_real_t  h_l)
{
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t t_l = h_l/(cs_glob_air_props->cp_l) - tkelvi;

  return  t_l;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the specific enthalpy of liquid water
 *
 * \return specific enthalpy of liquid water
 *
 * \param[in]     t_l            liquid water temperature (in Celsius)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_liq_t_to_h(cs_real_t  t_l)
{
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t h_l = cs_glob_air_props->cp_l * (t_l + tkelvi);

  return  h_l;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the air humidity at saturation for a given temperature
 *
 * \return absolute humidity of saturated air
 *
 * \param[in]     t_c            temperature in Celsius degree
 * \param[in]     p            reference pressure
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_x_sat(cs_real_t  t_c,
             cs_real_t  p)
{
  cs_real_t  pv;
  cs_real_t  x_s = 0.;

  /* Warning if T less than -20 degrees C */
  /* T between -20 and 80 degrees C */

  if ((t_c <= 80.)) {

    pv = cs_air_pwv_sat(t_c);
    x_s = 0.622 * pv/(p-pv);

  }

  /* T more than 80 degrees C */

  else if (t_c > 80.) {

    x_s = 0.5 + 0.001*t_c;

  }

  return x_s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the air water mass fraction at saturation for a
 * given temperature
 *
 * \return the air water mass fraction at saturation
 *
 * \param[in]  t_c   temperature in Celsius degree
 * \param[in]  p     reference pressure
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_yw_sat(cs_real_t  t_c,
              cs_real_t  p)
{
  cs_real_t  x_s = 0.;
  cs_real_t  y_s = 0.;

  x_s = cs_air_x_sat(t_c, p);
  y_s = x_s/(1. + x_s);

  return y_s;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the saturation water vapour pressure function of the
 *      temperature (in Celsius)
 *
 * \return the saturation water vapour pressure (=esatliq)
 *
 * \param[in]  t_c   temperature in Celsius degree
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_pwv_sat(cs_real_t  t_c)
{
  cs_real_t  a1, b1, c1, ps, pv;

  /* T between -20 and 0 degrees C */

  if (t_c <= 0.) {

    a1  = 6.4147;
    b1  = 22.376;
    c1  = 271.68;

    /* Warning if T less than -20 degrees C */

    ps  = a1 + (b1 * t_c)/(c1 + t_c);
    pv  = exp(ps);

  }

  /* T between 0 and 40 degrees C */

  else if (t_c <= 40.) {

    a1  = 6.4147;
    b1  = 17.438;
    c1  = 239.78;
    ps  = a1 + (b1 * t_c)/(c1 + t_c);
    pv  = exp(ps);
  }

  /* T between 40 and 80 degrees C */

  else {

    const cs_real_t  t0 = 273.16;
    const cs_real_t  ax = 8.2969;
    const cs_real_t  ay = 4.76955;
    const cs_real_t  a0 = 0.78614;
    const cs_real_t  a2 = 5.028;
    const cs_real_t  a3 = 0.000150475;
    const cs_real_t  a4 = 0.00042873;
    cs_real_t  tt, px, py, g1, g2, g3, g4;

    a1 = 10.7954;

    tt = t_c/t0;
    /* T greater than 80 degrees C, clipped at 80°C */
    if (t_c > 80.)
      tt = 80./t0;
    px = ax * tt;
    py = ay * tt/(1. + tt);
    g1 = a1 * tt/(1. + tt);
    g2 = -a2 * log10(1. + tt);
    g3 = a3 * (1. - 1./pow(10.,px));
    g4 = a4 * (pow(10., py) - 1.);
    ps = a0 + g1 + g2 + g3 + g4;
    pv = pow(10., ps) * 100.;

  }

  return pv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert the absolute humidity of humid air to the air
 *  water mass fraction qwt = Ym = mw/mh
 *
 * \return air water mass fraction
 *
 * \param[in]     x             absolute humidity of humid air
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_x_to_yw(cs_real_t  x)
{
  cs_real_t  qw = 0.;

  qw = x / (1. + x);

  return qw;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert the air water mass fraction qwt = Ym = mw/mh to the absolute
 * humidity of humid air
 *
 * \return absolute humidity of humid air
 *
 * \param[in]     qw             air water mass fraction
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_yw_to_x(cs_real_t  qw)
{
  cs_real_t  x = 0.;

  x = qw / (1. - qw);

  return x;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the density of humid air
 *
 * \param[in]     ywm           air water mass fraction
 * \param[in]     t_liq         temperature computed from
 *                              liquid potential temperature (K)
 * \param[in]     p             pressure
 * \param[out]    yw_liq        liquid water mass fraction
 * \param[out]    t_h           temperature of humid air in Celsius
 * \param[out]    rho_h         density of humid air
 */
/*----------------------------------------------------------------------------*/

void
cs_rho_humidair(cs_real_t   ywm, //TODO rename yw_h
                cs_real_t   t_liq,
                cs_real_t   p,
                cs_real_t  *yw_liq,
                cs_real_t  *t_h,
                cs_real_t  *rho_h)
{
  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  /* Mixture equivalent R perfect gas constant */
  cs_real_t lrhum;

  cs_real_t rair = phys_pro->r_pg_cnst;
  cs_real_t rvsra = phys_pro->rvsra;
  cs_real_t clatev = phys_pro->clatev;
  cs_real_t cp0 = phys_pro->cp0;

  /* Saturated vapor content */
  cs_real_t yw_sat = cs_air_yw_sat(*t_h, p);
  cs_real_t delta_yw = ywm - yw_sat;

  /* Temperature of the mixture in Kelvin */
  cs_real_t t_h_k = t_liq;

  /* Density of the air parcel
   * ------------------------- */

  /* Unsaturated air parcel */
  if (delta_yw <= 0.) {
    lrhum = rair*(1. + (rvsra - 1.)*ywm);
    *yw_liq = 0.;
  }
  /* Saturated (ie. with liquid water) air parcel */
  else {
    *yw_liq = delta_yw
      / (1. + yw_sat * cs_math_pow2(clatev)
      / (rair*rvsra*cp0*cs_math_pow2(t_h_k))); //FIXME rair * rvsra =rvap
    lrhum = rair*(1. + (rvsra - 1.)*(ywm - *yw_liq) -* yw_liq);
    t_h_k += (clatev/cp0) * *yw_liq;
  }

  /* Temperature of the mixture in Celsius */
  *t_h = t_h_k - cs_physical_constants_celsius_to_kelvin;
  /* Perfect gas law */
  *rho_h = p / (lrhum * t_h_k);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the density of humid air
 *
 * \return density of humid air
 *
 * \param[in]     x             absolute humidity of humid air
 * \param[in]     rho0          reference density of humid air
 * \param[in]     p0            reference pressure
 * \param[in]     t0            reference temperature of humid air
 * \param[in]     molmassrat    dry air to water vapour molecular mass ratio
 * \param[in]     t_h           temperature of humid air in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_rho_humidair(cs_real_t  x,
                    cs_real_t  rho0,
                    cs_real_t  p0,
                    cs_real_t  t0,
                    cs_real_t  molmassrat,
                    cs_real_t  t_h)
{
  cs_real_t  rho_h;

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t x_s = cs_air_x_sat(t_h, p0);

  if (x <= x_s)
    rho_h = rho0*(t0/(t_h+tkelvi))*molmassrat / (molmassrat+x);
  else {
    rho_h = rho0*(t0/(t_h+tkelvi))*molmassrat / (molmassrat+x_s);
    cs_real_t rho_l;
    if (t_h <= 0.)
      rho_l = 917.0;
    else
      rho_l = 998.36 - 0.4116 * (t_h-20.)
            - 2.24 * (t_h - 20.) * (t_h - 70.)/625.;//FIXME clipping if > 100 deg?

    rho_h  = 1. / (1. / rho_h + (x - x_s)/rho_l);
  }

  /* humid air */
  rho_h *= (1. + x);

  return rho_h;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
