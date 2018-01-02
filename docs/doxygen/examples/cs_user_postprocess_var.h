/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*!
  \page cs_user_postprocess_var Output additional variables on a post-processing mesh (cs_user_postprocess_var.f90)

  \section cs_user_postprocess_var_h_intro Introduction

  The \ref usvpst user subroutine allows one to output additional variables
  on a post-processing mesh. Several "automatic" post-processing meshes may be
  defined :
    - The volume mesh (\c ipart=-1)
    - The boundary mesh (\c ipart=-2)
    - SYRTHES coupling surface (\c ipart < -2)
    - Cooling tower exchange zone meshes (\c ipart < -2) if \c ichrze = 1

   Additional meshes (cells or faces) may also be defined through the GUI or
   using the \ref cs_user_postprocess_meshes function from the
   cs_user_postprocess.c file.

   The examples of post-processing given below are using the meshes defined
   \ref cs_user_postprocess "here".

  \section cs_user_postprocess_var_h_volume_mesh Output on the volume mesh (ipart = -1)

  \subsection cs_user_postprocess_var_h_volume_mesh_mom Output of a combination of moments

  A combination of moments can also be post-processed using the \ref usvpst subroutine.

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_2

  \section cs_user_postprocess_var_h_boundary_mesh Output on the boundary mesh (ipart = -2)

  Variables can be post-processed on the boundary mesh even if they are orignally
  define at cell centers. The following code block illustrates the post-processing
  of the density on the boundary mesh.

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_3

  \section cs_user_postprocess_var_h_user_mesh_1_2 Output on user meshes 1 or 2

  User meshes appearing in the examples are defined \ref cs_user_postprocess "here".

  \subsection cs_user_postprocess_var_h_user_mesh_1_2_vel Output of an interpolated velocity on user meshes

  An interpolated velocity is computed on both interior and boundary faces using a simple linear
  interpolation on user meshes 1 or 2.

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_4

  \subsection cs_user_postprocess_var_h_user_mesh_1_2_pr Output of the pressure on user meshes

  Similarly, the pressure is computed on both interior and boundary faces and then post-processed on user meshes
  1 or 2.

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_5

  \subsection cs_user_postprocess_var_h_user_mesh_cdg Output of the centers of gravity in different ways

  The examples below illustrate how to output a same variable in different
  ways (interlaced or not, using an indirection or not).

  \subsubsection cs_user_postprocess_var_h_user_mesh_cdg1 Output of the centers of gravity, interleaved

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_6

  \subsubsection cs_user_postprocess_var_h_user_mesh_cdg2 Output of the centers of gravity, non-interleaved, time-dependant

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_7

  \subsubsection cs_user_postprocess_var_h_user_mesh_cdg3 Output of the centers of gravity, with indirection (parent-based)

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_8

  \section cs_user_postprocess_var_h_user_mesh_3_4 Output on user meshes 3 or 4

  \subsection cs_user_postprocess_var_h_user_mesh_3_4_vel Output of an interpolated velocity on user meshes

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_9

  \subsection cs_user_postprocess_var_h_user_mesh_3_4_pr Output of the pressure on user meshes

  \snippet cs_user_postprocess_var.f90 postprocess_var_ex_10

*/
