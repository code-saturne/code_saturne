#ifndef __CS_RESTART_MAP_H__
#define __CS_RESTART_MAP_H__

/*============================================================================
 * Checkpoint / restart extension to mapped meshes
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Local type definitions
 *============================================================================*/


/*=============================================================================
 * Global variables
 *============================================================================*/


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate restart files should be mapped to a given mesh input
 *
 * \param[in]  mesh_path           path to mesh input
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_mesh_input(const char  *mesh_path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set options relative to restart file mapping to a given mesh input.
 *
 * \param[in]  tolerance_base      associated base tolerance (used for bounding
 *                                 box check only, not for location test)
 * \param[in]  tolerance_fraction  associated fraction of element bounding
 *                                 boxes added to tolerance
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_options(float  tolerance_base,
                           float  tolerance_fraction);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build mapping of restart files to different mesh if defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_build(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free restart file mapping to different mesh.
 *
 * Revert restart reading to default behavior.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_free(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RESTART_MAP_H__ */
