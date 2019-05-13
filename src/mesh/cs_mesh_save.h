#ifndef __CS_MESH_SAVE_H__
#define __CS_MESH_SAVE_H__

/*============================================================================
 * Save mesh Preprocessor data
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Save a mesh as preprocessor data.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   mb       <-- pointer to optional mesh builder structure, or NULL
 *   path     <-- optional directory name for output, or NULL for default
 *                (directory automatically created if necessary)
 *   filename <-- file name
 *----------------------------------------------------------------------------*/

void
cs_mesh_save(cs_mesh_t          *mesh,
             cs_mesh_builder_t  *mb,
             const char         *path,
             const char         *filename);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_SAVE_H__ */

