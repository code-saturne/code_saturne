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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_air_props.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Global variables
 *============================================================================*/

cs_ctwr_fluid_props_t  *cs_glob_ctwr_props = NULL;

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Definition des proprietes physiques
 *
 * Interface Fortran :
 *
 * SUBROUTINE CTPROF
 *----------------------------------------------------------------------------*/

void CS_PROCF (ctprof, CTPROF)
(
  const cs_real_t   *cpa,   /* Capacite calorifique de l air          */
  const cs_real_t   *cpv,   /* Capacite calorifique de la vapeur      */
  const cs_real_t   *cpe,   /* Capacite calorifique de l eau          */
  const cs_real_t   *hv0,   /* Chaleur latente                        */
  const cs_real_t   *rhoe,  /* Masse volumique de l eau               */
  const cs_real_t   *visc,  /* Viscosite Dynamique                    */
  const cs_real_t   *cond,  /* Conductivite                           */
  const cs_real_t   *gravx, /* gravite selon x                        */
  const cs_real_t   *gravy, /* gravite selon y                        */
  const cs_real_t   *gravz  /* gravite selon z                        */
)
{
  if (cs_glob_ctwr_props == NULL)
    BFT_MALLOC(cs_glob_ctwr_props, 1, cs_ctwr_fluid_props_t);

  cs_glob_ctwr_props->cpa   = *cpa;
  cs_glob_ctwr_props->cpv   = *cpv;
  cs_glob_ctwr_props->cpe   = *cpe;
  cs_glob_ctwr_props->hv0   = *hv0;
  cs_glob_ctwr_props->rhoe  = *rhoe;
  cs_glob_ctwr_props->visc  = *visc;
  cs_glob_ctwr_props->cond  = *cond;
  cs_glob_ctwr_props->gravx = *gravx;
  cs_glob_ctwr_props->gravy = *gravy;
  cs_glob_ctwr_props->gravz = *gravz;
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
 * DOUBLE PRECISION XSAT :  -> : absolute humidity of saturated air
 *----------------------------------------------------------------------------*/

void CS_PROCF (xsath, XSATH)
(
 const cs_real_t  *th,
       cs_real_t  *xsat
)
{
  *xsat = cs_ctwr_xsath(*th);
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
 * DOUBLE PRECISION DXSAT :  -> : derivative of the humidity of saturated air
 *----------------------------------------------------------------------------*/

void CS_PROCF (dxsath, DXSATH)
(
 const cs_real_t  *th,
       cs_real_t  *dxsat
)
{
  *dxsat =  cs_ctwr_dxsath(*th);
}


/*============================================================================
 * Public function definitions
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
cs_ctwr_xsath(const cs_real_t  th)
{
  const cs_real_t  Pres = 101325.;
  cs_real_t  a1,b1,c1,ps,pv;
  cs_real_t  xsat = 0.;

  /* T less than -20 degrees */

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
    xsat  = 0.622 * pv/(Pres - pv);

  }

  /* T between 0 and 40 degrees C */

  else if ((th >= 0.) && (th <= 40.)) {

    a1  = 6.4147;
    b1  = 17.438;
    c1  = 239.78;
    ps  = a1 + (b1 * th)/(c1 + th);
    pv  = exp(ps);
    xsat  = 0.622 * pv/(Pres - pv);
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
    xsat = 0.622 * pv/(Pres-pv);

  }

  else if (th > 80.) {

    xsat = 0.5 + 0.001*th;

  }

  return xsat;
}

/*----------------------------------------------------------------------------
 * Calculation of moist air mass enthalpy
 *
 * parameters:
 *   hair <-- absolute humidity of saturated air
 *   tair <-- air temperature in Celsius degree
 *
 * returns:
 *   air mass enthalpy
 *----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_enthair(const cs_real_t  xair,
                const cs_real_t  tair)
{
  const cs_real_t        Cpa = 1006. ;
  const cs_real_t        Cpv = 1831. ;
  const cs_real_t        Hvo = 2501600. ;
  cs_real_t  hair;

  hair  = ( Cpa +  xair * Cpv )* tair + xair * Hvo;

  return  hair;
}

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
cs_ctwr_heau(const cs_real_t  teau)
{
  const cs_real_t  Cpe =  4179.0 ;
  cs_real_t        heau;

  heau = Cpe * teau;

  return  heau;
}

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
cs_ctwr_dxsath(const cs_real_t  th)
{
  const cs_real_t  Pres = 101325.;
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

    dxsath = 0.622 * pv * Pres * grpim /pow(Pres - pv,2.);

  }

  /* T between 0 and 40 degrees C */

  else if (th >=0. && th <= 40.) {

    a1 = 6.4147;
    b1 = 17.438;
    c1 = 239.78;
    ps = a1 + (b1 * th)/(c1 + th);
    pv = exp(ps);
    grpim =  b1 *   c1 / pow(c1 + th,2.);

    dxsath =0.622 * pv * Pres * grpim /pow(Pres - pv,2.);

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
    dxsath= Pres * pvpr * 0.622 / pow (Pres - pv, 2);

  }
  else if (th > 80.) {

    dxsath =  0.001;

  }

  return dxsath ;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
