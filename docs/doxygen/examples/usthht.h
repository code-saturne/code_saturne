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
  \page us_thht Examples of enthalpy-temperature conversion law (usthht.f90)

   The \ref usthht subroutine is used to encapsulate a simple
   enthalpy-temperature conversion law and its inverse.

  \section usthht_h_intro Introduction

  The \ref usthht function allows one to define a simple
  enthalpy-temperature conversion law and its inverse. The
  parameters \c mode allows one to know in which direction
  the conversion will be made.

  \section usthht_h_usthht Examples

  The following code blocks show two examples of entahlpy-temperature
  conversion law.

  \subsection usthht_h_usthht_1 Example 1

  The first example corresponds to a simple law : \f[H = C_p  dt \f]

  \snippet usthht.f90 example_1

  \subsection usthht_h_usthht_2 Example 2

  The second example corresponds to a simple interpolation based on
  a tabulation defined hereafter and declared as a variable :

  \snippet usthht.f90 loc_var_dec

  The enthalpy-temperature conversion law is then defined :

  \snippet usthht.f90 example_2

*/
