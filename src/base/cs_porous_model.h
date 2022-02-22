#ifndef __CS_POROUS_MODEL_H__
#define __CS_POROUS_MODEL_H__

/*============================================================================
 * Porous model management
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

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/* Choice of the porous model */
extern int cs_glob_porous_model;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute fluid volumes and fluid surfaces in addition to cell volumes
 * and surfaces.
 *
 * parameters:
 *   porous_model <-- porous model option (> 0 for porosity)
 *----------------------------------------------------------------------------*/

void
cs_porous_model_set_model(int  porous_model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize disable_flag
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_init_disable_flag(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set (unset) has_disable_flag
 *
 * \param[in]  flag   1: on, 0: off
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_set_has_disable_flag(int  flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Init fluid quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_init_fluid_quantities(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic computation of the face porosity and factors.
 *
 * This is useful for the integral porous model.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_auto_face_porosity(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POROUS_MODEL_H__ */
