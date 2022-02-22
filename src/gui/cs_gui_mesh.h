#ifndef __CS_GUI_MESH_H__
#define __CS_GUI_MESH_H__

/*============================================================================
 * Management of the GUI parameters file: mesh related options
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Determine whether warped faces should be cut.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_warping(void);

/*-----------------------------------------------------------------------------
 * Define joinings using a GUI-produced XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_define_joinings(void);

/*-----------------------------------------------------------------------------
 * Define periodicities using a GUI-produced XML file.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_define_periodicities(void);

/*----------------------------------------------------------------------------
 * Mesh smoothing.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure to smoothe
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_smoothe(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Define user thin wall through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_boundary(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Define user mesh extrude through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_extrude(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Define mesh save behavior trough the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_save_if_modified(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Define if cartesian mesh is to be built through GUI.
 *----------------------------------------------------------------------------*/

int
cs_gui_mesh_build_cartesian(void);

/*----------------------------------------------------------------------------
 * Read cartesian mesh parameters defined with GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_cartesian_define(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_MESH_H__ */
