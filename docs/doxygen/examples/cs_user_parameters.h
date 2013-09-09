/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
  \page parameters Input of calculation parameters (Fortran modules)
 
  \section intro Introduction

  User subroutines for input of calculation parameters (Fortran modules).
    These subroutines are called in all cases.

  If the Code_Saturne GUI is used, this file is not required (but may be
    used to override parameters entered through the GUI, and to set
    parameters not accessible through the GUI).

  Several routines are present in the file, each destined to defined
    specific parameters.

  To modify the default value of parameters which do not appear in the
    examples provided, code should be placed as follows:
    - usipsu   for numerical and physical options
    - usipes   for input-output related options

  As a convention, "specific physics" defers to the following modules only:
    pulverized coal, gas combustion, electric arcs.

  In addition, specific routines are provided for the definition of some
    "specific physics" options.
    These routines are described at the end of this file and will be activated
    when the corresponding option is selected in the usppmo routine.


  \section usipes  Input-output related options

  Frequency of log output.

  \snippet cs_user_parameters-output.f90 init_01

  Log (listing) verbosity.

  \snippet cs_user_parameters-output.f90 init_02

  Probes output step.

  \snippet cs_user_parameters-output.f90 init_03

  Number of monitoring points (probes) and their positions.
  Limited to ncaptm=100.

  \snippet cs_user_parameters-output.f90 init_04

  Current variable.
  As for other variables, if we do not assign the following array values,
  default values will be used:

     - ichrvr( ) = chonological output (yes 1/no 0)
     - ilisvr( ) = logging in listing (yes 1/no 0)
     - ihisvr( ) = history output (number of probes and their numbers)
     - if ihisvr(.,1)  = -1, output for all probes

  \note Only the fist 8 characters of a name will be used in the most detailed log.

  \snippet cs_user_parameters-output.f90 init_05

  User scalar variables.

  We may modify here the arrays relative to user scalars, but scalars
    reserved for specific physics are handled automatically. This explains
    the tests on 'nscaus', which ensure that the targeted scalars are
    truly user scalars.
  By specific physics, we mean only those which are handled in specific
    modules of the code, such as coal, combustion, electric arcs (see usppmo).

  \snippet cs_user_parameters-output.f90 init_06

  Other variables.

  \snippet cs_user_parameters-output.f90 init_07

  Specific physics variables.

  \snippet cs_user_parameters-output.f90 init_08

  Variables for coal particles.

  \snippet cs_user_parameters-output.f90 init_09

  Variables for droplets.

  \snippet cs_user_parameters-output.f90 init_10

  Variables for carrying gas.

  \snippet cs_user_parameters-output.f90 init_11

  Variables of State; User defined Variables.

  \snippet cs_user_parameters-output.f90 init_12

  State variables for coal particles or fuel droplets.

  \snippet cs_user_parameters-output.f90 init_13

  State variables for carrier gas phase.

  \snippet cs_user_parameters-output.f90 init_14


  \section examples Examples

  \subsection example_1 Example 1

   Force postprocessing of projection of some variables at boundary
   with no reconstruction.
   This is handled automatically if the second bit of a field's
   'post_vis' key value is set to 1 (which amounts to adding 2
   to that key value).

   field_get_id returns -1 if field does not exist

  
  \snippet cs_user_parameters-output.f90 example_1
  
  \subsection example_2 Example 2

  Enforce existence of 'tplus' and 'tstar' fields, so that
  a boundary temperature or Nusselt number may be computed using the
  post_boundary_temperature or post_boundary_nusselt subroutines.
  When postprocessing of these quantities is activated, those fields
  are present, but if we need to compute them in the
  cs_user_extra_operations user subroutine without postprocessing them,
  forcing the definition of these fields to save the values computed
  for the boundary layer is necessary.

  
  \snippet cs_user_parameters-output.f90 example_2
  
*/
