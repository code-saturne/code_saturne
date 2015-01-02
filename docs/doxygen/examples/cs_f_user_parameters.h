/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
    - usipsu   for numerical and physical options
    - usipes   for input-output related options

  As a convention, "specific physics" defers to the following modules only:
    pulverized coal, gas combustion, electric arcs.

  In addition, specific routines are provided for the definition of some
    "specific physics" options.
    These routines are described at the end of this file and will be activated
    when the corresponding option is selected in the usppmo routine.

  \section cs_f_user_parameters_h_usipsu  General options (usipsu)

  \subsection cs_f_user_parameters_h_usipsu_1 Special fields

  Enforce existence of 'tplus' and 'tstar' fields, so that
  a boundary temperature or Nusselt number may be computed using the
  post_boundary_temperature or post_boundary_nusselt subroutines.
  When postprocessing of these quantities is activated, those fields
  are present, but if we need to compute them in the
  \ref cs_user_extra_operations user subroutine without postprocessing them,
  forcing the definition of these fields to save the values computed
  for the boundary layer is necessary.

  \snippet cs_user_parameters-output.f90 usipsu_ex_1

  Save contribution of slope test for variables in special fields.
  These fields are automatically created, with postprocessing output enabled,
  if the matching variable is convected, does not use a pure upwind scheme,
  and has a slope test (the slope_test_upwind_id key value for a given
  variable's field is automatically set to the matching postprocessing field's
  id, or -1 if not applicable).

  \snippet cs_user_parameters-output.f90 usipsu_ex_2

  \section cs_f_user_parameters_h_usipes  Input-output related examples (usipes)

  \subsection cs_f_user_parameters_h_example_base Basic options

  Frequency of log output.

  \snippet cs_user_parameters-output.f90 usipes_ex_01

  Log (listing) verbosity.

  \snippet cs_user_parameters-output.f90 usipes_ex_02

  Activate or deactivate logging output.
  By default, logging is active for most variables. In the following
  example, logging for velocity is deactivated.

  \snippet cs_user_parameters-output.f90 usipes_ex_03

  Change a property's label (here for density, first checking if it
  is variable). A field's name cannot be changed, but its label,
  used for logging and postprocessing output, may be redefined.

  \snippet cs_user_parameters-output.f90 usipes_ex_04

  \subsection cs_f_user_parameters_h_example_probes Probes output

  Probes output step.

  \snippet cs_user_parameters-output.f90 usipes_ex_05

  Number of monitoring points (probes) and their positions.
  Limited to ncaptm=100.

  \snippet cs_user_parameters-output.f90 usipes_ex_06

  \subsection cs_f_user_parameters_h_example_post Postprocessing output

  Activate or deactivate postprocessing output.
  By default, output is active for most variables. In the following
  example, the output for velocity is deactivated.

  \snippet cs_user_parameters-output.f90 usipes_ex_07

  Activate or deactivate probes output.
  If \ref ihisvr(.,1) = -1, output is done for all probes.
  In the following example, probes output for the velocity is restricted
  to the first component.

  \snippet cs_user_parameters-output.f90 usipes_ex_08


	Probes for Radiative Transfer (Luminance and radiative density flux vector)
  for all probes (ihisvr = -1)

  \snippet cs_user_parameters-output.f90 usipes_ex_10

  \subsection cs_f_user_parameters_h_example_1 Postprocess at boundary

   Force postprocessing of projection of some variables at boundary
   with no reconstruction.
   This is handled automatically if the second bit of a field's
   'post_vis' key value is set to 1 (which amounts to adding 2
   to that key value).

  \snippet cs_user_parameters-output.f90 usipes_ex_09

*/
