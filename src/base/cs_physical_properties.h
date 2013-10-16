#ifndef __CS_PHYSICAL_PROPERTIES_H__
#define __CS_PHYSICAL_PROPERTIES_H__

/*============================================================================
 * Functions dealing with parallelism
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

typedef enum {

  PLANE_PH,
  PLANE_PT,
  PLANE_PS,
  PLANE_PU,
  PLANE_PV,
  PLANE_TS,
  PLANE_TX,

} cs_thermo_plane_type_t;

typedef enum {

  PRESSURE,
  TEMPERATURE,
  ENTHALPY,
  ENTROPY,
  ISOBARIC_HEAT_CAPACITY,
  ISOCHORIC_HEAT_CAPACITY,
  SPECIFIC_VOLUME,
  DENSITY,
  INTERNAL_ENERGY,
  QUALITY,
  THERMAL_CONDUCTIVITY,
  DYNAMIC_VISCOSITY,
  SPEED_OF_SOUND

} cs_property_type_t;


/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute properties with freeteam in a defined thermal plane.
 *
 * Fortran Interface :
 *
 * subroutine  csfspp(thermo_plane, property, ncel, var1, var2, val)
 * *****************
 *
 * integer           thermo_plane : <-- : define thermal plane
 * integer           property     : <-- : define property to compute
 * integer           ncel         : <-- : size of the local array
 * double precision  var1(*)      : <-- : first array for thermal plane
 * double precision  var2(*)      : <-- : second array for thermal plane
 * double precision  val(*)       : --> : values
 *----------------------------------------------------------------------------*/

void
CS_PROCF (csfspp, CSFSPP)(cs_thermo_plane_type_t *thermo_plane,
                          cs_property_type_t *property,
                          const    int *const ncel,
                          double *var1,
                          double *var2,
                          double *val);


END_C_DECLS

#endif /* __CS_PHYSICAL_PROPERTIES_H__ */
