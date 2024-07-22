#ifndef __CS_AIR_PROPS_H__
#define __CS_AIR_PROPS_H__

/*============================================================================
 * Specific laws for air properties (temperature, enthalpy)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Structure definition
 *============================================================================*/

/* Structure associated to general properties */

typedef struct {

  cs_real_t   humidity0;       /* Reference humidity */
  cs_real_t   cp_a;            /* Specific heat of dry air */
  cs_real_t   cp_v;            /* Specific heat of vapor */
  cs_real_t   cp_l;            /* Specific heat of liquid water */
  cs_real_t   hv0;             /* Chaleur latente */
  cs_real_t   rho_l;           /* Masse volumique de l eau */
  cs_real_t   lambda_h;        /* Humid air conductivity */
  cs_real_t   lambda_l;        /* Water conductivity */
  cs_real_t   droplet_diam;    /* Drop diameter for rain zones */
  cs_real_t   molmass_rat;     /* Ratio of the molar mass (H2O) over the
                                  molar mass (air) */
  cs_real_t   sigma;           /* Surface tension between water and air */

} cs_air_fluid_props_t;

extern  cs_air_fluid_props_t  *cs_glob_air_props;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*============================================================================
 *  Prototypes of public function
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
cs_air_h_l(cs_real_t  t_l);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation water vapor mass enthalpy
 *
 * \return water vapor mass enthalpy
 *
 * \param[in]     t_vap          water vapor temperature in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_hvap(cs_real_t  t_vap);

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
              cs_real_t  p0);

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
                   cs_real_t  x_s);

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
                  cs_real_t  t_h);

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
                  cs_real_t  h_h);

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
cs_liq_h_to_t(cs_real_t  h_l);

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
cs_liq_t_to_h(cs_real_t  t_l);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the saturation water vapor pressure function of the
 *      temperature (in Celsius)
 *
 * \return the saturation water vapor pressure (=esatliq)
 *
 * \param[in]  t_c   temperature in Celsius degree
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
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
 * \brief Calculation of the air humidity at saturation for a given temperature
 *
 * \return absolute humidity of saturated air
 *
 * \param[in]     t_c          temperature in Celsius degree
 * \param[in]     p            reference pressure
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
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
              cs_real_t  p);

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
cs_air_x_to_yw(const cs_real_t  x);

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
cs_air_yw_to_x(cs_real_t  qw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the density of humid air
 *
 * \param[in]     yw_h          air water mass fraction
 * \param[in]     theta_l       liquid potential temperature (K)
 * \param[in]     p             pressure
 * \param[out]    yw_liq        liquid water mass fraction
 * \param[out]    t_h           temperature of humid air in Celsius
 * \param[out]    rho_h         density of humid air
 * \param[out]    beta_h        thermal expansion of the bulk
 */
/*----------------------------------------------------------------------------*/

void
cs_rho_humidair(cs_real_t   yw_h,
                cs_real_t   t_liq,
                cs_real_t   p,
                cs_real_t  *yw_liq,
                cs_real_t  *t_h,
                cs_real_t  *rho_h,
                cs_real_t  *beta_h);

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
 * \param[in]     molmassrat    dry air to water vapor molecular mass ratio
 * \param[in]     t_h           temperature of humid air in Celsius
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_rho_humidair(cs_real_t  x,
                    cs_real_t  rho0,
                    cs_real_t  p0,
                    cs_real_t  t0,
                    cs_real_t  molmassrat,
                    cs_real_t  t_h);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_AIR_PROPERTIES_H__ */
