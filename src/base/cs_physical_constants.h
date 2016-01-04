#ifndef __CS_PHYSICAL_CONSTANTS_H__
#define __CS_PHYSICAL_CONSTANTS_H__

/*============================================================================
 * Base physical constants data.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* physical constants descriptor */
/*-------------------------------*/

typedef struct {

  double        gx, gy, gz;                 /* gravity components */
  int           icorio;                     /* Coriolis source terms indicator */

} cs_physical_constants_t;

/* fluid properties descriptor */
/*-----------------------------*/

typedef struct {

  int           ixyzp0;       /* indicator for filling of reference point for
                                 total pressure */
  double        ro0;          /* reference density */
  double        viscl0;       /* reference molecular dynamic viscosity */
  double        p0;           /* reference pressure for the total pressure */
  double        pred0;        /* reference value for the reduced pressure */
  double        xyzp0[3];     /* reference point coordinates for the total
                                 pressure */
  double        t0;           /* reference temperature */
  double        cp0;          /* reference specific heat */
  double        xmasmr;       /* molar mass of the perfect gas in kg/mol
                                 (if ieos=1) */
  double        pther;        /* uniform thermodynamic pressure for the low-Mach
                                 algorithm */
  double        pthera;       /* thermodynamic pressure for the previous time
                                 step */
  double        pthermax;     /* thermodynamic maximum pressure for user
                                 clipping, used to model a venting effect */

} cs_fluid_properties_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main physical constants structure */

extern const cs_physical_constants_t  *cs_glob_physical_constants;

/* Pointer to main fluid properties structure */

extern const cs_fluid_properties_t  *cs_glob_fluid_properties;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PHYSICAL_CONSTANTS_H__ */
