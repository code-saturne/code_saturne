/*============================================================================
 * Specific laws for air properties (temperature, enthalpy)
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
#include "cs_physical_constants.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_air_props.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main physical constants structure and associated pointer */

static cs_ctwr_fluid_props_t _props = {
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

cs_ctwr_fluid_props_t *cs_glob_ctwr_props = &_props;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

void
cs_ctwr_glob_properties_get_pointer(double  **humidity0,
                                    double  **cp_a,
                                    double  **cp_v,
                                    double  **cp_l,
                                    double  **hv0,
                                    double  **rho_l,
                                    double  **lambda_h,
                                    double  **lambda_l,
                                    double  **droplet_diam)
{

  *humidity0    = &(_props.humidity0    );
  *cp_a         = &(_props.cp_a         );
  *cp_v         = &(_props.cp_v         );
  *cp_l         = &(_props.cp_l         );
  *hv0          = &(_props.hv0          );
  *rho_l        = &(_props.rho_l        );
  *lambda_h     = &(_props.lambda_h     );
  *lambda_l     = &(_props.lambda_l     );
  *droplet_diam = &(_props.droplet_diam );
}

/*----------------------------------------------------------------------------
 * Calculation of the air humidity at saturation for a given temperature
 *
 * Fortran interface:
 *
 * SUBROUTINE XSATH
 * ****************
 *
 * DOUBLE PRECISION TH   : <-  : temperature in Celsius degree
 * DOUBLE PRECISION P0   : <-  : reference pressure
 * DOUBLE PRECISION XSAT :  -> : absolute humidity of saturated air
 *----------------------------------------------------------------------------*/

void CS_PROCF (xsath, XSATH)
(
 const cs_real_t  *th,
 const cs_real_t  *p0,
       cs_real_t  *xsat
)
{
  *xsat = cs_ctwr_xsath(*th, *p0);
}

/*----------------------------------------------------------------------------
 * Calculation of the derivative of the absolute humidity at saturation
 *
 * Fortran interface:
 *
 * SUBROUTINE DXSATH
 * *****************
 *
 * DOUBLE PRECISION TH    : <-  : temperature in Celsius degree
 * DOUBLE PRECISION P0    : <-  : reference pressure
 * DOUBLE PRECISION DXSAT :  -> : derivative of the humidity of saturated air
 *----------------------------------------------------------------------------*/

void CS_PROCF (dxsath, DXSATH)
(
 const cs_real_t  *th,
 const cs_real_t  *p0,
       cs_real_t  *dxsat
)
{
  *dxsat =  cs_ctwr_dxsath(*th, *p0);
}

/*----------------------------------------------------------------------------
 * Calculation of the density of humid air
 *
 * Fortran interface:
 *
 * SUBROUTINE RHO_HUMIDAIR
 * **********************
 *
 *   DOUBLE PRECISION HUMID        : <-  : absolute humidity of humid air
 *   DOUBLE PRECISION R0           : <-  : reference density of humid air
 *   DOUBLE PRECISION T0           : <-  : reference temperature of humid air
 *   DOUBLE PRECISION MOLMASSRAT   : <-  : dry air to water vapour molecular mass ratio
 *   DOUBLE PRECISION T_H          : <-  : temperature of humid air in Celsius
 *
 *   DOUBLE PRECISION RHO_HUMIDAIR : ->  : density of humid air
 *----------------------------------------------------------------------------*/

void CS_PROCF (rho_humidair, RHO_HUMIDAIR)
(
 const cs_real_t *x,
 const cs_real_t *rho0,
 const cs_real_t *p0,
 const cs_real_t *t0,
 const cs_real_t *molmassrat,
 const cs_real_t *t_h,
       cs_real_t *rho_humidair
)
{
  *rho_humidair =  cs_ctwr_rho_humidair(*x, *rho0, *p0, *t0, *molmassrat, *t_h);
}

/*----------------------------------------------------------------------------
 * Calculation of the Cp of humid air
 *
 * Fortran interface:
 *
 * SUBROUTINE CP_HUMIDAIR
 * **********************
 *
 *   DOUBLE PRECISION HUMID       : <-  : absolute humidity of humid air
 *   DOUBLE PRECISION HUMID_SAT   : <-  : absolute humidity of saturated humid air
 *
 *   DOUBLE PRECISION CP_HUMIDAIR : ->  : specific heat of humid air
 *----------------------------------------------------------------------------*/

void CS_PROCF (cp_humidair, CP_HUMIDAIR)
(
 const cs_real_t  *x,
 const cs_real_t  *x_s,
       cs_real_t  *cp_humidair
)
{
  *cp_humidair =  cs_ctwr_cp_humidair(*x,*x_s);
}

/*----------------------------------------------------------------------------
 * Calculation of the specific enthalpy of humid air
 *
 * Fortran interface:
 *
 * SUBROUTINE H_HUMIDAIR
 * *********************
 *
 *   DOUBLE PRECISION CP_HUMIDAIR : <-  : specific heat of humid air
 *   DOUBLE PRECISION X           : <-  : absolute humidity of humid air
 *   DOUBLE PRECISION HUMID_SAT   : <-  : absolute humidity of saturated humid air
 *   DOUBLE PRECISION T_HUMIDAIR  : <-  : humid air temperature in Celsius (in Celsius)
 *
 *   DOUBLE PRECISION H_HUMIDAIR  : ->  : specific enthalpy of humid air
 *----------------------------------------------------------------------------*/

void CS_PROCF (h_humidair, H_HUMIDAIR)
(
 const cs_real_t  *cp_humidair,
 const cs_real_t  *x,
 const cs_real_t  *x_s,
 const cs_real_t  *t_humidair,
       cs_real_t  *h_humidair
)
{
  *h_humidair =  cs_ctwr_h_humidair(*cp_humidair,*x,*x_s,*t_humidair);
}

/*----------------------------------------------------------------------------
 * Calculation of the temperature of humid air
 *
 * Fortran interface:
 *
 * SUBROUTINE T_HUMIDAIR
 * *********************
 *
 *   DOUBLE PRECISION CP_HUMIDAIR : <-  : specific heat of humid air
 *   DOUBLE PRECISION X           : <-  : absolute humidity of humid air
 *   DOUBLE PRECISION HUMID_SAT   : <-  : absolute humidity of saturated humid air
 *   DOUBLE PRECISION H_HUMIDAIR  : <-  : specific enthalpy of humid air
 *
 *   DOUBLE PRECISION T_HUMIDAIR  : ->  : humid air temperature in Celsius (in Celsius)
 *----------------------------------------------------------------------------*/

void CS_PROCF (t_humidair, T_HUMIDAIR)
(
 const cs_real_t  *cp_humidair,
 const cs_real_t  *x,
 const cs_real_t  *x_s,
 const cs_real_t  *h_humidair,
       cs_real_t  *t_humidair

)
{
  *t_humidair =  cs_ctwr_t_humidair(*cp_humidair,*x,*x_s,*h_humidair);
}

/*----------------------------------------------------------------------------
 * Calculation of the specific enthalpy of liquid water
 *
 * Fortran interface:
 *
 * SUBROUTINE H_LIQUIDWATER
 * ************************
 *
 *   DOUBLE PRECISION T_LIQWATER : <-  : liquid water temperature (in Celsius)
 *
 *   DOUBLE PRECISION H_LIQWATER : ->  : specific enthalpy of liquid water
 *----------------------------------------------------------------------------*/

void CS_PROCF (h_liqwater, H_LIQWATER)
(
 const cs_real_t  *t_liqwater,
       cs_real_t  *h_liqwater
)
{
  *h_liqwater =  cs_ctwr_h_liqwater(*t_liqwater);
}

/*----------------------------------------------------------------------------
 * Calculation of the specific enthalpy of humid air at a
 * specified humidity and temperature
 *
 * Fortran interface:
 *
 * SUBROUTINE H_HUMIDAIR_FIXED
 * ***************************
 *
 *   DOUBLE PRECISION T_AIR : <-  : humid air temperature (in Celsius)
 *   DOUBLE PRECISION X_AIR : <-  : humid air humidity
 *
 *   DOUBLE PRECISION H_AIR : ->  : specific enthalpy of humid air
 *----------------------------------------------------------------------------*/

void CS_PROCF (h_humidair_fixed, H_HUMIDAIR_FIXED)
(
 const cs_real_t  *x_air,
 const cs_real_t  *t_air,
 cs_real_t  *h_humidair_fixed
)
{
  *h_humidair_fixed =  cs_ctwr_enthair(*x_air,*t_air);//FIXME rename because it is generic
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the air humidity at saturation for a given temperature
 *
 * \return absolute humidity of saturated air
 *
 * \param[in]     th            temperature in Celsius degree
 * \param[in]     p0            reference pressure
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_xsath(const cs_real_t  th,
              const cs_real_t  p0)
{
  cs_real_t  a1,b1,c1,ps,pv;
  cs_real_t  xsat = 0.;

  /* T less than -20 degrees C */

  if (th  < -20.) {

    xsat = 0.;

  }

  /* T between -20 and 0 degrees C */

  else if ((th >= -20.) && (th <= 0.)) {

    a1  = 6.4147;
    b1  = 22.376;
    c1  = 271.68;
    ps  = a1 + (b1 * th)/(c1 + th);
    pv  = exp(ps);
    xsat  = 0.622 * pv/(p0 - pv);

  }

  /* T between 0 and 40 degrees C */

  else if ((th >= 0.) && (th <= 40.)) {

    a1  = 6.4147;
    b1  = 17.438;
    c1  = 239.78;
    ps  = a1 + (b1 * th)/(c1 + th);
    pv  = exp(ps);
    xsat  = 0.622 * pv/(p0 - pv);
  }

  /* T between 40 and 80 degrees C */

  else if ((th >= 40.) && (th <= 80.)) {

    const cs_real_t  T0 = 273.16;
    const cs_real_t  Ax = 8.2969;
    const cs_real_t  Ay = 4.76955;
    const cs_real_t  A0 = 0.78614;
    const cs_real_t  A2 = 5.028;
    const cs_real_t  A3 = 0.000150475;
    const cs_real_t  A4 = 0.00042873;
    cs_real_t  tt,px,py,g1,g2,g3,g4;

    a1 = 10.7954;

    tt = th/T0;
    px = Ax * tt;
    py = Ay * tt/(1. + tt);
    g1 = a1 * tt/(1. + tt);
    g2 = -A2 * log10(1. + tt);
    g3 = A3 * (1. - 1./pow(10.,px));
    g4 = A4 * (pow(10.,py) - 1.);
    ps = A0 + g1 + g2 + g3 + g4;
    pv = pow(10.,ps) * 100.;
    xsat = 0.622 * pv/(p0-pv);

  }

  else if (th > 80.) {

    xsat = 0.5 + 0.001*th;

  }

  return xsat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of moist air mass enthalpy
 *
 * \return air mass enthalpy
 *
 * \param[in]     xair          absolute humidity of saturated air
 * \param[in]     tair          air temperature in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_enthair(const cs_real_t  xair,
                const cs_real_t  tair)
{
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  cs_real_t  h_air;

  h_air  = (ct_prop->cp_a + xair*ct_prop->cp_v) * tair
         + xair*ct_prop->hv0;//FIXME /1+x a priori not

  return  h_air;
}

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
cs_ctwr_h_l(const cs_real_t  t_l)
{
  cs_real_t h_l = (cs_glob_ctwr_props->cp_l) * t_l;

  return  h_l;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation water vapour mass enthalpy
 *
 * \return water vapour mass enthalpy
 *
 * \param[in]     t_vap          water vapour temperature in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_hvap(const cs_real_t  t_vap)
{
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

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
cs_ctwr_dxsath(const cs_real_t  th,
               const cs_real_t  p0)
{
  cs_real_t   a1,b1,c1,ps,pv,grpim;
  cs_real_t   dxsath = 0.;

  /* T less than -20 degrees */

  if (th < -20.) {

    dxsath = 0.;

  }

  /* T between -20 and 0 degrees C */

  else if (th >= -20. && th <= 0.) {

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
cs_ctwr_rho_humidair(const cs_real_t x,
                     const cs_real_t rho0,
                     const cs_real_t p0,
                     const cs_real_t t0,
                     const cs_real_t molmassrat,
                     const cs_real_t t_h)
{
  cs_real_t  rho_h;

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t x_s = cs_ctwr_xsath(t_h, p0);

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
cs_ctwr_cp_humidair(const cs_real_t  x,
                    const cs_real_t  x_s)
{
  cs_real_t  cp_h;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  if (x <= x_s) {
    cp_h = ct_prop->cp_a + x*ct_prop->cp_v;
  } else {
    cp_h = ct_prop->cp_a + x_s*ct_prop->cp_v + (x-x_s)*ct_prop->cp_l;
  } //FIXME + dX_s * hv ?

  cp_h = cp_h/(1.0+x);

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
cs_ctwr_h_humidair(const cs_real_t  cp_h,
                   const cs_real_t  x,
                   const cs_real_t  x_s,
                   const cs_real_t  t_h)
{
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;
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
cs_ctwr_t_humidair(const cs_real_t  cp_h,
                   const cs_real_t  x,
                   const cs_real_t  x_s,
                   const cs_real_t  h_h)
{
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;
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
cs_ctwr_t_liqwater(const cs_real_t h_l)
{
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t t_l = h_l/(cs_glob_ctwr_props->cp_l) - tkelvi;

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
cs_ctwr_h_liqwater(const cs_real_t t_l)
{
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t h_l = cs_glob_ctwr_props->cp_l * (t_l + tkelvi);

  return  h_l;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
