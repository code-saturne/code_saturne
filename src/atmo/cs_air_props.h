#ifndef __CS_AIR_PROPS_H__
#define __CS_AIR_PROPS_H__

/*============================================================================
 * Specific laws for air properties (temperature, enthalpy)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
  cs_real_t   cp_a;            /* Capacite calorifique de l air */
  cs_real_t   cp_v;            /* Capacite calorifique de la vapeur */
  cs_real_t   cp_l;            /* Capacite calorifique de l eau */
  cs_real_t   hv0;             /* Chaleur latente */
  cs_real_t   rho_l;           /* Masse volumique de l eau */
  cs_real_t   lambda_h;        /* Humid air conductivity */
  cs_real_t   lambda_l;        /* Water conductivity */
  cs_real_t   droplet_diam;    /* Drop diameter for rain zones */

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
 * \return water vapour mass enthalpy
 *
 * \param[in]     t_vap          water vapour temperature in Celsius
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
             cs_real_t  p);

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
 * \brief Computes the saturation water vapour pressure function of the
 *      temperature (in Celsius)
 *
 * \return the saturation water vapour pressure (=esatliq)
 *
 * \param[in]  t_c   temperature in Celsius degree
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_air_pwv_sat(cs_real_t  t_c);

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
                cs_real_t  *rho_h);

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
                    cs_real_t  t_h);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_AIR_PROPERTIES_H__ */
