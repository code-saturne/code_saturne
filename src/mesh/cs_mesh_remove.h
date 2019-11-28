#ifndef __CS_MESH_REMOVE_H__
#define __CS_MESH_REMOVE_H__

/*============================================================================
 * Functions to remove mesh elements.
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

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove flagged cells.
 *
 * \param[in, out]  m            mesh
 * \param[in]       flag        cell flag (!= 0 to remove)
 * \param[in]       group_name  name of group to assign to new boundary faces,
 *                              or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_remove_cells(cs_mesh_t    *m,
                     char          flag[],
                     const char   *group_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Remove cells with negative volumes
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_remove_cells_negative_volume(cs_mesh_t  *m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_REMOVE_H__ */
