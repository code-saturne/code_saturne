#ifndef __CS_LAGR_ROUGHNESS_H__
#define __CS_LAGR_ROUGHNESS_H__

/*============================================================================
 * Functions and types for the roughness modeling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  cs_real_t   faraday_cst;
  cs_real_t   free_space_permit;
  cs_real_t   water_permit;
  cs_real_t   ionic_strength;
  cs_real_t   phi1;
  cs_real_t   phi2;
  cs_real_t  *temperature;
  cs_real_t  valen;
  cs_real_t  *debye_length;
  cs_real_t   cstham;
  cs_real_t   dcutof;
  cs_real_t   lambwl;
  cs_real_t   kboltz;
  cs_real_t   espasg;
  cs_real_t   denasp;
  cs_real_t   rayasp;
  cs_real_t   rayasg;

} cs_lagr_roughness_param_t;

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Clogging initialization:
 *  - Retrieve various parameters for storing in global structure
 *  - Compute and store the Debye screening length
 *----------------------------------------------------------------------------*/

void
CS_PROCF (roughness_init, ROUGHNESS_INIT)(const cs_real_t   *faraday_cst,
                                          const cs_real_t   *free_space_permit,
                                          const cs_real_t   *water_permit,
                                          const cs_real_t   *ionic_strength,
                                          const cs_real_t    temperature[],
                                          const cs_real_t   *valen,
                                          const cs_real_t   *phi1,
                                          const cs_real_t   *phi2,
                                          const cs_real_t   *cstham,
                                          const cs_real_t   *dcutof,
                                          const cs_real_t   *lambwl,
                                          const cs_real_t   *kboltz,
                                          const cs_real_t   *espasg,
                                          const cs_real_t   *denasp,
                                          const cs_real_t   *rayasp,
                                          const cs_real_t   *rayasg
);

/*=============================================================================
 * Function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Deallocate the arrays storing temperature and Debye length.
 *----------------------------------------------------------------------------*/

void
cs_lagr_roughness_finalize(void);

/*----------------------------------------------------------------------------
 * Compute the energy barrier for a rough wall.
 *
 * parameters:
 *   particle       <-- pointer to particle data
 *   attr_map       <-- pointer to attribute map
 *   face_id        <-- id of face neighboring particle
 *   energy_barrier <-> energy barrier
 *----------------------------------------------------------------------------*/

void
cs_lagr_roughness_barrier(const void                     *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lnum_t                       face_id,
                          cs_real_t                      *energy_barrier);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_ROUGHNESS_H__ */

