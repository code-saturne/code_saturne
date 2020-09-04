#ifndef __CS_STL_H__
#define __CS_STL_H__

/*============================================================================
 * STL reader and writer.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an STL mesh
 *----------------------------------------------------------------------------*/

typedef struct {

  char          name[10];      /*!< Name identifier of the STL file*/

  char          header[80];    /*!< Header of the STL file */

  cs_lnum_t     n_faces;       /*!< Number of triangles */

  float        *coords;        /*!< Array of face vertex coordinates
                                *  vtx_coord[n_vertices][3] */

  fvm_nodal_t  *ext_mesh;      /*!< Associated external mesh */

} cs_stl_mesh_t ;

/*----------------------------------------------------------------------------
 * Structure containing the number of STL meshes and their associated pointers
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_stl_mesh_t  **mesh_list;  /*!< Array of STL meshes
                                *   size: n_meshes*/

  cs_lnum_t        n_meshes;   /*!< Total number of STL meshes */

} cs_stl_mesh_info_t ;

/*=============================================================================
 * Static global variables
 *============================================================================*/

extern cs_stl_mesh_info_t  *cs_glob_stl_meshes;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Add a new STL mesh structure.
 *
 * parameters:
 *   path <-- path of the STL mesh
 *
 * returns:
 *   pointer to the new STL mesh structure
 *----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_add(const char  *path);

/*----------------------------------------------------------------------------
 * Return a pointer to a STL mesh based on its name if present.
 *
 * parameters:
 *   name <-- name of the STL mesh
 *
 * returns:
 *   pointer to the STL mesh structure, or NULL
 *----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_get_by_name(const char  *name);

/*----------------------------------------------------------------------------
 * Free all allocated STL mesh structures
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_destroy_all(void);

/*----------------------------------------------------------------------------
 * Read a binary STL file and store its content in a STL mesh structure.
 *
 * parameters:
 *   stl_mesh  <-- pointer to the associated STL mesh structure
 *   path      <-- path to the STL file
 *   matrix    <-- transformation matrix
 *----------------------------------------------------------------------------*/

void
cs_stl_file_read(cs_stl_mesh_t  *stl_mesh,
                 const char     *path,
                 double         matrix[3][4]);

/*----------------------------------------------------------------------------
 * Write a binary STL file according to a given STL mesh structure.
 *
 * parameters:
 *   stl_mesh  <-- pointer to the associated STL mesh structure
 *   path      <-- path to the STL file
 *----------------------------------------------------------------------------*/

void
cs_stl_file_write(cs_stl_mesh_t  *stl_mesh,
                  const char     *path);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_STL_H__ */
