/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_LAGR_CLOGGING_H__
#define __CS_LAGR_CLOGGING_H__

/*============================================================================
 * Functions and types for the clogging modeling
 *============================================================================*/

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
  cs_real_t   jamming_limit;
  cs_real_t   min_porosity;
  cs_real_t   phi1;
  cs_real_t   phi2;

  cs_real_t*  temperature;
  cs_real_t*  debye_length;
  cs_real_t   cstham;
  cs_real_t   dcutof;
  cs_real_t   lambwl;
  cs_real_t   kboltz

} cs_lagr_clog_param_t;


/*----------------------------------------------------------------------------
 * Clogging initialization:
 *  - Retrieve various parameters for storing in global structure
 *  - Compute and store the Debye screening length
 *----------------------------------------------------------------------------*/

void
CS_PROCF (cloginit, CLOGINIT)(const cs_real_t   *faraday_cst,
                              const cs_real_t   *free_space_permit,
                              const cs_real_t   *water_permit,
                              const cs_real_t   *ionic_strength,
                              const cs_real_t   *jamming_limit,
                              const cs_real_t   *min_porosity,
                              const cs_real_t    temperature[],
                              const cs_real_t   *phi1,
                              const cs_real_t   *phi2,
                              const cs_real_t   *cstham,
                              const cs_real_t   *dcutof,
                              const cs_real_t   *lambwl,
                              const cs_real_t   *kboltz
);


/*----------------------------------------------------------------------------
 * Clogging ending
 * Deallocate the arrays storing temperature and Debye length
 *----------------------------------------------------------------------------*/

void clogend( void );


/*----------------------------------------------------------------------------
 * Clogging:
 *
 * - Calculation of the number of deposited particles in contact with the
 *   depositing particle
 * - Re-calculation of the energy barrier in this number is greater than zero
 *----------------------------------------------------------------------------*/


cs_int_t
clogging_barrier(cs_lagr_particle_t     particle,
                 cs_int_t               face_id,
                 cs_real_t              face_area,
                 cs_real_t*             energy_barrier,
                 cs_real_t*             surface_coverage,
                 cs_real_t*             limit ,
                 cs_real_t*             mporos

  );

cs_real_t
vdwsa(           cs_real_t              distcc,
                 cs_real_t              rpart1,
                 cs_real_t              rpart2
  );


cs_real_t
edlsa(            cs_int_t               ifac,
                  cs_real_t              distcc,
                  cs_real_t              rpart1,
                  cs_real_t              rpart2
                 );


END_C_DECLS

#endif /* __CS_LAGR_CLOGGING_H__ */

