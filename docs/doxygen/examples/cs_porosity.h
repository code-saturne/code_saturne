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
  \page cs_porosity Examples of data settings for porous media (cs_user_porosity.f90)


  \section intro_poro Introduction

  This function computes the porosity (volume factor \f$ \epsilon \f$
  when porosity module is activated (iporos = 1 in cs_user_parameters.f90).

  \section cs_user_poro_examples Porosity setting examples
  Here is the list of examples:

  - \subpage base_poro_examples

*/
// __________________________________________________________________________________
/*!

  \page base_poro_examples Basic example


  \section base_loc_var_poro Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_porosity.f90 loc_var_dec


  \section subroutine_end_poro Test to remove for use of subroutine end

  The following initialization block needs to be added for the following examples:


  \subsection porosity_field Retrieve porosity field

  \snippet cs_user_porosity.f90 init  


  \subsection example_porosity Example: fixe a linear by part porosity profile

  \snippet cs_user_porosity.f90 example_1

  Periodicity and parallelism treatment:

  \snippet cs_user_porosity.f90 parallelism


*/
// __________________________________________________________________________________
