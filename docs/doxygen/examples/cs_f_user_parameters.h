/*============================================================================
 * Code_Saturne documentation page
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

/*-----------------------------------------------------------------------------*/

/*!
  \page f_parameters Input of calculation parameters (Fortran modules)

  \section cs_f_user_parameters_h_intro Introduction

  User subroutines for input of calculation parameters (Fortran modules).
    These subroutines are called in all cases.

  If the Code_Saturne GUI is used, this file is not required (but may be
    used to override parameters entered through the GUI, and to set
    parameters not accessible through the GUI).

  Several routines are present in the file, each destined to defined
    specific parameters.

  To modify the default value of parameters which do not appear in the
    examples provided, code should be placed as follows:
    - \ref cs_f_user_parameters_h_usipsu "usipsu" for numerical and physical options
    - \ref cs_f_user_parameters_h_usipes "usipes" for input-output related options
    - \ref cs_f_user_parameters_h_usppmo "usppmo" for specific physics options
    - \ref cs_f_user_parameters_h_usipph "usipph" for additional input of parameters
    - \ref cs_f_user_parameters_h_usati1 "usati1" for calculation options for the atmospheric module
    - \ref cs_f_user_parameters_h_cs_user_combustion "cs_user_combustion" for calculation options for the combustion module
    - \ref cs_f_user_parameters_h_uscfx1 "uscfx1" and \ref cs_f_user_parameters_h_uscfx2 "uscfx2" for non-standard options for the compressible module
    - \ref cs_f_user_parameters_h_uscti1 "uscti1" for the definition of cooling tower model and exchange zones
    - \ref cs_f_user_parameters_h_user_darcy_ini1 "user_darcy_ini1" for calculation options for the Darcy module

  As a convention, "specific physics" defers to the following modules only:
    pulverized coal, gas combustion, electric arcs.

  In addition, specific routines are provided for the definition of some
    "specific physics" options.
    These routines are described at the end of this file and will be activated
    when the corresponding option is selected in the \ref usppmo routine.

  \section cs_f_user_parameters_h_usipsu  General options (usipsu)

  \subsection cs_f_user_parameters_h_usipsu_0 All options

  The following code block presents all the options available
  in the \ref usipsu subroutine.

  \snippet cs_user_parameters.f90 usipsu

  \section cs_f_user_parameters_h_usppmo Specific physic activation (usppmo)

  The \ref usppmo routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 usppmo

  \section cs_f_user_parameters_h_usipph Additional input of parameters (usipph)

  The \ref usipph routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 usipph

  \section cs_f_user_parameters_h_usati1 Calculation options for the atmospheric module (usati1)

  The \ref usati1 routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 usati1

  \section cs_f_user_parameters_h_cs_user_combustion Calculation options for the combustion module (cs_user_combustion)

  The \ref cs_user_combustion routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 cs_user_combustion

  \section cs_f_user_parameters_h_uscfx1 Non-standard options for the compressible module (uscfx1)

  The \ref uscfx1 routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 uscfx1

  \section cs_f_user_parameters_h_uscfx2 Non-standard options for the compressible module (uscfx2)

  The \ref uscfx2 routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 uscfx2

  \section cs_f_user_parameters_h_uscti1 Definition of cooling tower model and exchange zones (uscti1)

  The \ref uscti1 routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 uscti1

  \section cs_f_user_parameters_h_user_darcy_ini1 Calculation options for the Darcy module (user_darcy_ini1)

  The \ref user_darcy_ini1 routine can be found in the \ref cs_user_parameters.f90 file.

  \snippet cs_user_parameters.f90 user_darcy_ini1

*/
