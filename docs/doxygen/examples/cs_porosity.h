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

/*-----------------------------------------------------------------------------*/

/*!
  \page cs_porosity Examples of data settings for porous media
  (cs_user_porosity.c)

  \section intro_poro Introduction

  This function computes the porosity (volume factor \f$ \epsilon \f$
  when porosity module is activated (iporos = 1 in cs_user_parameters.f90).

  \section cs_user_poro_examples Porosity setting examples
  Here is the list of examples:

  - \subpage base_poro_examples

*/

/*-----------------------------------------------------------------------------*/

/*!
  \page base_poro_examples Setting porosity values: basic example

  \section base_loc_var_poro Local definitions and initialization

  \subsection mesh_quantities Mesh quantities

  It may be useful to access some mesh adjacencies and  quantities,
  in which case local pointers allow for more readable code:

  \snippet cs_user_porosity.c init_poro_mq

  \subsection properties Associated properties

  Accessing cell porosity property values is required so values may be set:

  \snippet cs_user_porosity.c init_poro_pro

  \section example_porosity Example: define porosity by geometric zones

  Individual cell porosity values can be assigned to each cell, so
  they may be based on groups, geometric criteria, or any other
  time-independent functions:

  \snippet cs_user_porosity.c set_poro_cells_1

  Matching face equivalent surfaces should also be assigned in a
  corresponding manner, for interior faces:

  \snippet cs_user_porosity.c set_poro_i_faces_1

  and for boundary faces:

  \snippet cs_user_porosity.c set_poro_b_faces_1
*/

/*-----------------------------------------------------------------------------*/
