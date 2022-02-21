#ifndef __CS_WALLCONDENSATION_1DTHERMAL_H__
#define __CS_WALLCONDENSATION_1DTHERMAL_H__

/*============================================================================
 * Base wall condensation model data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

typedef struct {
  int        nzones;
  int        znmurx;
  cs_real_t *ztheta;
  cs_real_t *zdxmin;
  cs_lnum_t *znmur;
  cs_real_t *zepais;
  cs_real_t *ztpar0;

  cs_real_t *zhext;
  cs_real_t *ztext;
  cs_real_t *zrob;
  cs_real_t *zcondb;
  cs_real_t *zcpb;
  cs_real_t *ztpar;

  cs_real_t *zdxp;
  cs_real_t *ztmur;

} cs_wall_cond_1d_thermal_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to wall condensation descriptor structure */

extern const cs_wall_cond_1d_thermal_t *cs_glob_wall_cond_1d_thermal;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for wall condensation thermal models.
 *
 * \param[in] nzones number of zones
 */
/*----------------------------------------------------------------------------*/

void cs_wall_condensation_1d_thermal_create(int nzones);

void
cs_wall_condensation_1d_thermal_mesh_create(int znmurx, int nfbpcd, int nzones);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to wall condensation models.
 */
/*----------------------------------------------------------------------------*/

void cs_wall_condensation_1d_thermal_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide writeable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_cond_1d_thermal_t *cs_get_glob_wall_cond_1d_thermal(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALLCONDENSATION_1DTHERMAL_H__ */
