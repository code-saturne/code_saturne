#ifndef __CS_LAGR_DLVO_H__
#define __CS_LAGR_DLVO_H__

/*============================================================================
 * Functions and types for the clogging modeling
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  cs_real_t   water_permit;
  cs_real_t   ionic_strength;
  cs_real_t   phi_p;
  cs_real_t   phi_s;
  cs_real_t  *temperature;
  cs_real_t   valen;
  cs_real_t  *debye_length;
  cs_real_t   cstham;
  cs_real_t   csthpp;
  cs_real_t   lambda_vdw;

} cs_lagr_dlvo_param_t;

/*=============================================================================
 * Function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * DLVO initialization
 *----------------------------------------------------------------------------*/

void
cs_lagr_dlvo_init(const cs_real_t   water_permit,
                  const cs_real_t   ionic_strength,
                  const cs_real_t   temperature[],
                  const cs_real_t   valen,
                  const cs_real_t   phi_p,
                  const cs_real_t   phi_s,
                  const cs_real_t   cstham,
                  const cs_real_t   csthpp,
                  const cs_real_t   lambda_vdw);

/*----------------------------------------------------------------------------
 * Deallocate the arrays storing temperature and Debye length.
 *----------------------------------------------------------------------------*/

void
cs_lagr_dlvo_finalize(void);

/*----------------------------------------------------------------------------
 * Compute the energy barrier for a smooth wall.
 *----------------------------------------------------------------------------*/

void
cs_lagr_barrier(const void                     *particle,
                const cs_lagr_attribute_map_t  *attr_map,
                cs_lnum_t                       iel,
                cs_real_t                      *energy_barrier);

/*----------------------------------------------------------------------------
 * Compute the energy barrier for two spheres.
 *----------------------------------------------------------------------------*/

void
cs_lagr_barrier_pp(cs_real_t                       dpart,
                   cs_lnum_t                       iel,
                   cs_real_t                      *energy_barrier);

/*----------------------------------------------------------------------------
 * Van der Waals interaction between a sphere and a plane
 * using formulas from Czarnecki (large distances)
 *                 and Gregory (small distances)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_plane(cs_real_t  distp,
                                   cs_real_t  rpart,
                                   cs_real_t  lambda_vdw,
                                   cs_real_t  cstham);

/*----------------------------------------------------------------------------
 * Calculation of the Van der Waals interaction between two spheres
 * following the formula from Gregory (1981a)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_sphere(cs_real_t    distcc,
                                    cs_real_t    rpart1,
                                    cs_real_t    rpart2,
                                    cs_real_t    lambda_vdw,
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
                          cs_real_t  temp,
                          cs_real_t  debye_length,
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
                          cs_real_t  temp,
                          cs_real_t  debye_length,
                          cs_real_t  water_permit);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_DLVO_H__ */

