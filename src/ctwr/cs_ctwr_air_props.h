#ifndef __CS_CTWR_AIR_PROPS_H__
#define __CS_CTWR_AIR_PROPS_H__

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

} cs_ctwr_fluid_props_t;

extern  cs_ctwr_fluid_props_t  *cs_glob_ctwr_props;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Calculation of the air humidity at saturation for a given temperature
 *
 * Fortran interface:
 *
 * SUBROUTINE XSATH
 * ****************
 *
 * DOUBLE PRECISION TH   : <-  : temperature in Celsius degree
 * DOUBLE PRECISION XSAT :  -> : absolute humidity of saturated air
 *----------------------------------------------------------------------------*/

void CS_PROCF (xsath, XSATH)
(
 const cs_real_t  *th,
 const cs_real_t  *p0,
       cs_real_t  *xsat
);

/*----------------------------------------------------------------------------
 * Calculation of the derivative of the absolute humidity at saturation
 *
 * Fortran interface:
 *
 * SUBROUTINE DXSATH
 * *****************
 *
 * DOUBLE PRECISION TH    : <-  : temperature in Celsius degree
 * DOUBLE PRECISION DXSAT :  -> : derivative of the humidity of saturated air
 *----------------------------------------------------------------------------*/

void CS_PROCF (dxsath, DXSATH)
(
 const cs_real_t  *th,
 const cs_real_t  *p0,
       cs_real_t  *dxsat
);

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
 const cs_real_t *humid,
 const cs_real_t *r0,
 const cs_real_t *p0,
 const cs_real_t *t0,
 const cs_real_t *molmassrat,
 const cs_real_t *t_h,
       cs_real_t *rho_humidair
);

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
 const cs_real_t  *humid,
 const cs_real_t  *humid_sat,
       cs_real_t  *cp_humidair
);

/*----------------------------------------------------------------------------
 * Calculation of the specific enthalpy of humid air
 *
 * Fortran interface:
 *
 * SUBROUTINE H_HUMIDAIR
 * *********************
 *
 *   DOUBLE PRECISION CP_HUMIDAIR : <-  : specific heat of humid air
 *   DOUBLE PRECISION HUMID       : <-  : absolute humidity of humid air
 *   DOUBLE PRECISION HUMID_SAT   : <-  : absolute humidity of saturated humid air
 *   DOUBLE PRECISION T_HUMIDAIR  : <-  : humid air temperature in Celsius (in Celsius)
 *
 *   DOUBLE PRECISION H_HUMIDAIR  : ->  : specific enthalpy of humid air
 *----------------------------------------------------------------------------*/

void CS_PROCF (h_humidair, H_HUMIDAIR)
(
 const cs_real_t  *cp_humidair,
 const cs_real_t  *humid,
 const cs_real_t  *humid_sat,
 const cs_real_t  *t_humidair,
       cs_real_t  *h_humidair
);

/*----------------------------------------------------------------------------
 * Calculation of the temperature of humid air
 *
 * Fortran interface:
 *
 * SUBROUTINE T_HUMIDAIR
 * *********************
 *
 *   DOUBLE PRECISION CP_HUMIDAIR : <-  : specific heat of humid air
 *   DOUBLE PRECISION HUMID       : <-  : absolute humidity of humid air
 *   DOUBLE PRECISION HUMID_SAT   : <-  : absolute humidity of saturated humid air
 *   DOUBLE PRECISION H_HUMIDAIR  : <-  : specific enthalpy of humid air
 *
 *   DOUBLE PRECISION T_HUMIDAIR  : ->  : humid air temperature in Celsius (in Celsius)
 *----------------------------------------------------------------------------*/

void CS_PROCF (t_humidair, T_HUMIDAIR)
(
 const cs_real_t  *cp_humidair,
 const cs_real_t  *humid,
 const cs_real_t  *humid_sat,
 const cs_real_t  *h_humidair,
       cs_real_t  *t_humidair
);

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

void CS_PROCF (h_liqwater, H_LIQUWATER)
(
 const cs_real_t  *t_liqwater,
       cs_real_t  *h_liqwater
);

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

void CS_PROCF (h_humidair_fixed, H_HUMIDAIR_FIXED)(
  const cs_real_t  *x_air,
  const cs_real_t  *t_air,
  cs_real_t        *h_humidair_fixed);

/*============================================================================
 *  Prototypes of public function
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Calculation of the air humidity at saturation for a given temperature
 *
 * parameters:
 *   th <-- temperature in Celsius degree
 *
 * returns:
 *   absolute humidity of saturated air
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_xsath(const cs_real_t  th,
              const cs_real_t  p0);

/*----------------------------------------------------------------------------
 * Calculation of moist air mass enthalpy
 *
 * parameters:
 *   xair <-- absolute humidity of saturated air
 *   tair <-- air temperature in Celsius degree
 *
 * returns:
 *   air mass enthalpy
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_enthair(const cs_real_t  xair,
                const cs_real_t  tair);

/*----------------------------------------------------------------------------
 * Calculation water mass enthalpy
 *
 * parameters:
 *   teau <-- water temperature in Celsius degree
 *
 * returns:
 *   water mass enthalpy
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_heau(const cs_real_t  teau);

/*----------------------------------------------------------------------------
 * Calculation water vapour mass enthalpy
 *
 * parameters:
 *   tvap <-- water vapour temperature in Celsius degree
 *
 * returns:
 *   water vapour mass enthalpy
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_hvap(const cs_real_t  t_vap);

/*----------------------------------------------------------------------------
 * Calculation of the derivate of the absolute humidity at saturation
 *
 * parameters:
 *   th <-- temperature in Celsius degree
 *
 * returns:
 *   derivative of the humidity of saturated air
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_dxsath(const cs_real_t  th,
               const cs_real_t  p0);

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
                     const cs_real_t t_h);

/*----------------------------------------------------------------------------
 * Calculation of the Cp of humid air
 *
 * parameters:
 *   humid     <-- absolute humidity of humid air
 *   humid_sat <-- absolute humidity of saturated humid air
 *
 * returns:
 *   cp_h <-- specific heat of humid air
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_cp_humidair(const cs_real_t  humid,
                    const cs_real_t  humid_sat);

/*----------------------------------------------------------------------------
 * Calculation of the specific enthalpy of humid air
 *
 * parameters:
 *   cp_h      <-- Cp of humid air
 *   humid     <-- absolute humidity of humid air
 *   humid_sat <-- absolute humidity of saturated humid air
 *   t_h       <-- humid air temperature in Celsius (in Celsius)
 *
 * returns:
 *   h_h <-- specific enthalpy of humid air
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_h_humidair(const cs_real_t  cp_h,
                   const cs_real_t  humid,
                   const cs_real_t  humid_sat,
                   const cs_real_t  t_h);

/*----------------------------------------------------------------------------
 * Calculation of the temperature of humid air
 *
 * parameters:
 *   cp_h      <-- Cp of humid air
 *   humid     <-- absolute humidity of humid air
 *   humid_sat <-- absolute humidity of saturated humid air
 *   h_h       <-- humid air enthalpy
 *
 * returns:
 *   t_h <-- temperature of humid air (in Celsius)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_t_humidair(const cs_real_t  cp_h,
                   const cs_real_t  humid,
                   const cs_real_t  humid_sat,
                   const cs_real_t  h_h);

/*----------------------------------------------------------------------------
 * Calculation of the temperature of liquid water
 *
 * parameters:
 *   h_liqwater <-- specific enthalpy of liquid water
 *
 * returns:
 *   liquid water temperature (in Celsius)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_t_liqwater(const cs_real_t h_liqwater);

/*----------------------------------------------------------------------------
 * Calculation of the specific enthalpy of liquid water
 *
 * parameters:
 *   t_l <-- liquid water temperature (in Celsius)
 *
 * returns:
 *   specific enthalpy of liquid water
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_h_liqwater(const cs_real_t t_liqwater);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_AIR_PROPERTIES_H__ */
