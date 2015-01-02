#ifndef __CS_LAGR_DLVO_H__
#define __CS_LAGR_DLVO_H__

/*============================================================================
 * Functions and types for the clogging modeling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Van der Waals interaction between a sphere and a plane
 * using formulas from Czarnecki (large distances)
 *                 and Gregory (small distances)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_plane(cs_real_t  distp,
                                   cs_real_t  rpart,
                                   cs_real_t  lambwl,
                                   cs_real_t  cstham);

/*----------------------------------------------------------------------------
 * Calculation of the Van der Waals interaction between two spheres
 * following the formula from Gregory (1981a)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_sphere(cs_real_t    distcc,
                                    cs_real_t    rpart1,
                                    cs_real_t    rpart2,
                                    cs_real_t    lambwl,
                                    cs_real_t    cstham);

/*----------------------------------------------------------------------------
 * Electric Double Layer (EDL) interaction between a sphere and a plane
 * using the formula from Bell & al (1970)
 * based on the McCartney & Levine method
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_edl_sphere_plane (cs_real_t  distp,
                          cs_real_t  rpart,
                          cs_real_t  valen,
                          cs_real_t  phi1,
                          cs_real_t  phi2,
                          cs_real_t  kboltz,
                          cs_real_t  temp,
                          cs_real_t  debye_length,
                          cs_real_t  free_space_permit,
                          cs_real_t  water_permit);

/*----------------------------------------------------------------------------
 * Calculation of the EDL interaction between two spheres
 * using the formula from Bell & al (1970)
 * based on the McCartney & Levine method
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_edl_sphere_sphere(cs_real_t  distcc,
                          cs_real_t  rpart1,
                          cs_real_t  rpart2,
                          cs_real_t  valen,
                          cs_real_t  phi1,
                          cs_real_t  phi2,
                          cs_real_t  kboltz,
                          cs_real_t  temp,
                          cs_real_t  debye_length,
                          cs_real_t  free_space_permit,
                          cs_real_t  water_permit);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_DLVO_H__ */

