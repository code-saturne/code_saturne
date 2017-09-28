#ifndef __CS_PHYSICAL_CONSTANTS_H__
#define __CS_PHYSICAL_CONSTANTS_H__

/*============================================================================
 * Base physical constants data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

  cs_real_3_t   gravity;           /* gravity vector */
  int           icorio;            /* Coriolis source terms indicator */

} cs_physical_constants_t;

/* fluid properties descriptor */
/*-----------------------------*/

typedef struct {

  int           ixyzp0;       /* indicator for filling of reference point for
                                 total pressure */
  int           ieos;         /* indicator of equation of state for the
                                 compressible module */
  int           icp;          /* property index of the isobaric specific heat */
  int           icv;          /* property index of the isochoric specific
                                 heat */
  int           irovar;       /* variable density field */
  int           ivivar;       /* variable viscosity field */
  int           ivsuth;       /* Sutherland law for laminar viscosity and
                                 thermal conductivity in gas mix spec. phys. */
  double        ro0;          /* reference density */
  double        viscl0;       /* reference molecular dynamic viscosity */
  double        p0;           /* reference pressure for the total pressure */
  double        pred0;        /* reference value for the reduced pressure */
  double        xyzp0[3];     /* reference point coordinates for the total
                                 pressure */
  double        t0;           /* reference temperature */
  double        cp0;          /* reference specific heat */
  double        cv0;          /* reference specific heat */
  double        xmasmr;       /* molar mass of the perfect gas in kg/mol
                                 (if ieos=1) */
  double        psginf;       /* stiffened gas (ieos=2) limit pressure (zero in
                                 perfect gas) (Pa) */
  double        gammasg;      /* stiffened gas (ieos=2) polytropic coefficient,
                                 (dimensionless) */
  int           ipthrm;       /* Uniform variable thermodynamic pressure for the
                                 low-Mach algorithm */
  double        pther;        /* uniform thermodynamic pressure for the low-Mach
                                 algorithm */
  double        pthera;       /* thermodynamic pressure for the previous time
                                 step */
  double        pthermax;     /* thermodynamic maximum pressure for user
                                 clipping, used to model a venting effect */
  double        sleak;        /* Leak surface */
  double        kleak;        /* Leak head loss (2.9 by default, from Idelcick) */
  double        roref;        /* Initial reference density */

} cs_fluid_properties_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Physical constants */

extern const double cs_physical_constants_r; /* Ideal gas constant (J/mol/K) */
extern const double cs_physical_constants_kb; /* Boltzmann constant (J/K) */
extern const double cs_physical_constants_celsius_to_kelvin; /* Celsius to
                                                                Kelvin*/
extern const double cs_physical_constants_stephan; /* Stephan constant
                                                     (W/m2/K4)*/

/* Pointer to main physical constants structure */

extern const cs_physical_constants_t  *cs_glob_physical_constants;

/* Pointer to main fluid properties structure */

extern const cs_fluid_properties_t  *cs_glob_fluid_properties;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_physical_constants
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_physical_constants_t *
cs_get_glob_physical_constants(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_fluid_properties
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_fluid_properties_t *
cs_get_glob_fluid_properties(void);

/*----------------------------------------------------------------------------
 * Print the physical constants structure to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_physical_constants_log_setup(void);

/*----------------------------------------------------------------------------
 * Print the fluid properties structure to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_fluid_properties_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PHYSICAL_CONSTANTS_H__ */
