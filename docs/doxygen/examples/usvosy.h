/*============================================================================
 * Code_Saturne documentation page
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

/*-----------------------------------------------------------------------------*/

/*!
  \page us_vosy Examples of volume exchange coefficient computation for SYRTHES coupling (usvosy.f90)

   The \ref usvosy subroutine is used to compute a volume exchange
   coefficient for SYRTHES coupling.

  \section usvosy_h_usvosy Examples

  The following code blocks show two examples of computation of
  a volume exchange coefficient.

  \subsection usvosy_h_usvosy_arg Arguments of usvosy

  \snippet usvosy.f90 arg

  \subsection usvosy_h_usvosy_loc_var_dec Variable declaration

  \snippet usvosy.f90 loc_var_dec

  \subsection usvosy_h_usvosy_init Initialization

  The values of the different fields that will be needed for the computation of the
  volume exchange coefficient are retrieved.

  \snippet usvosy.f90 init

  \subsection usvosy_h_usvosy_1 Example 1

  The first example corresponds to a constant volume exchange coefficient.

  \snippet usvosy.f90 example_1

  \subsection usvosy_h_usvosy_2 Example 2

  The second example corresponds to a variable volume exchange coefficient
  defined as follows :

  \f[ h_{vol} = h_{surf} S \f]

  with S is the surface area where exchanges take place by unit of volume and

  \f[ h_{surf} = \frac{Nu \lambda}{L} \f]

  \snippet usvosy.f90 example_2

*/
