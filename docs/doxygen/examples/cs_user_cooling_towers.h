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
  \page cs_user_cooling_towers Cooling tower parameters definition

  \section cs_user_cooling_towers_h_fluid Define fluid parameters

  Fluid properties are defined in \ref cs_user_parameters.f90 in routine
  \ref cs_user_cooling_towers:

  \snippet cs_user_parameters.f90 cs_user_cooling_towers

  \section cs_user_cooling_towers_h_ex Define an evaporation zone

  The \ref cs_user_paramters.c function allows one to define parameters
  for cooling towers.

  The following code block shows an example of definition of a cooling tower
  exchange zone.

  \snippet cs_user_parameters-ctwr.c ctwr_user_1

*/
