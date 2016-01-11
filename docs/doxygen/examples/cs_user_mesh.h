/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*-----------------------------------------------------------------------------*/

/*!
  \page cs_user_mesh Mesh modification

  \section cs_user_mesh_h_intro Introduction

  C user functions for optional modification of the mesh.
    These subroutines are called in all cases.

  Several functions are present in the file, each specific to different
    modification types.

  \section cs_user_mesh_h_cs_user_mesh_modifiy  General mesh modifications

  Mesh modifications not available through specialized functions
  should be defined in \ref cs_user_mesh_modify.

  \subsection cs_user_mesh_h_cs_user_mesh_modifiy_coords Coordinates modification

  For example, to modify coordinates, the
  following code can be added:

  \snippet cs_user_mesh-modify.c mesh_modify_coords

  \subsection cs_user_mesh_h_cs_user_mesh_modifiy_extrude_1 Boundary mesh patch extrusion

  It is possible to add cells by extruding selected boundary faces.
  A simplified usage is available through the \ref cs_mesh_extrude_constant,
  while extruding a mesh with vertex-local extrusion parameters is available
  through \ref cs_mesh_extrude.

  The example which follows illustrate using the simplified function.

  \snippet cs_user_mesh-modify.c mesh_modify_extrude_1

*/
