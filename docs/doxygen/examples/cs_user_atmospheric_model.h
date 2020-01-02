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
  \page cs_user_atmospheric_model Atmospheric model (cs_user_atmospheric_model.f90)

  \section cs_user_atmospheric_model_h_intro Introduction

  User subroutine for the atmospheric model.

  \section cs_user_atmospheric_model_usatdv Atmospheric module

  \c imode corresponds to the number of calls of the \ref usatdv function.
  Depending on the value of \c imode , different operations are performed
  in the following example.

  \subsection cs_user_atmospheric_model_usatdv_imode0 imode = 0

  \snippet cs_user_atmospheric_model.f90 imode_0

  \subsection cs_user_atmospheric_model_usatdv_imode1 imode = 1

  \snippet cs_user_atmospheric_model.f90 imode_1

  \section cs_user_atmospheric_model_usatsoil Data Entry for the atmospheric ground model

  \snippet cs_user_atmospheric_model.f90 usatsoil

*/
