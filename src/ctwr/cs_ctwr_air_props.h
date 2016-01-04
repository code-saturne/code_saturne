#ifndef __CS_CTWR_AIR_PROPS_H__
#define __CS_CTWR_AIR_PROPS_H__

/*============================================================================
 * Specific laws for air properties (temperature, enthalpy)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

  cs_real_t  cpa;             /* Capacite calorifique de l air */
  cs_real_t  cpv;             /* Capacite calorifique de la vapeur */
  cs_real_t  cpe;             /* Capacite calorifique de l eau */
  cs_real_t  hv0;             /* Chaleur latente */
  cs_real_t  rhoe;            /* Masse volumique de l eau*/
  cs_real_t  visc;            /* Viscosite Dynamique */
  cs_real_t  cond;            /* Conductivite */
  cs_real_t  gravx;           /* Gravite x */
  cs_real_t  gravy;           /* Gravite y */
  cs_real_t  gravz;           /* Gravite z */

} cs_ctwr_fluid_props_t;

extern  cs_ctwr_fluid_props_t  *cs_glob_ctwr_props;

/* Structure associated to air properties */

typedef struct {

  cs_real_t  rho_ref; /* Reference density */
  cs_real_t  p_ref;   /* Reference pressure */
  cs_real_t  t_ref;   /* Reference temperature */
  cs_real_t  delta;

  cs_real_t  g[3]; /* Gravity vector */

} cs_ctwr_air_props_t;

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
       cs_real_t  *dxsat
);

/*----------------------------------------------------------------------------
 * Communication des proprietes physiques
 *
 * Fortran interface:
 *
 * SUBROUTINE COMPPF
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (ctprof, CTPROF)
(
 const cs_real_t  *cpa,             /* Capacite calorifique de l air */
 const cs_real_t  *cpv,            /* Capacite calorifique de la vapeur */
 const cs_real_t  *cpe,            /* Capacite calorifique de l eau */
 const cs_real_t  *hv0,            /* Chaleur latente */
 const cs_real_t  *rhoe,           /* Masse volumique de l eau*/
 const cs_real_t  *visc,           /* Viscosite Dynamique */
 const cs_real_t  *cond,           /* Conductivite */
 const cs_real_t  *gravx,          /* Gravite x */
 const cs_real_t  *gravy,          /* Gravite y */
 const cs_real_t  *gravz           /* Gravite z */
);

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
cs_ctwr_xsath(const cs_real_t  th);

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
 * Calculation of the derivate of the absolute humidity at saturation
 *
 * parameters:
 *   th <-- temperature in Celsius degree
 *
 * returns:
 *   derivative of the humidity of saturated air
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_dxsath(const cs_real_t  th);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_AIR_PROPERTIES_H__ */
