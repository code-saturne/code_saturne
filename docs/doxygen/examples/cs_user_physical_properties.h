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
  \page physical_properties Define user laws for physical properties

  \section cs_user_physical_properties_h_intro Introduction

  This page provides examples of code blocks that may be used
  to define physical variable laws.

  \warning

  It is \b forbidden to modify turbulent viscosity \c visct here
  (a specific subroutine is dedicated to that: \ref usvist)

  - icp = 1 must <b> have been specified </b>
     in \ref usipsu if we wish to define a variable specific heat
     cpro_cp (otherwise: memory overwrite).

  - the kivisl field integer key (diffusivity_id)
    must <b> have been specified </b>
    in \ref usipsu if we wish to define a variable viscosity
    \c viscls.

  \remarks
   - This routine is called at the beginning of each time step
     Thus, <b> AT THE FIRST TIME STEP </b> (non-restart case), the only
     values initialized before this call are those defined
       - in the GUI or  \ref usipsu (cs_user_parameters.f90)
              - density    (initialized at \c ro0)
              - viscosity  (initialized at \c viscl0)
       - in the GUI or \ref cs_user_initialization
              - calculation variables (initialized at 0 by defaut
              or to the value given in the GUI or in \ref cs_user_initialization)

   - We may define here variation laws for cell properties, for:
      - density:                                    rom    kg/m3
      - density at boundary faces:                  romb   kg/m3)
      - molecular viscosity:                        cpro_viscl  kg/(m s)
      - specific heat:                              cpro_cp     J/(kg degrees)
      - diffusivities associated with scalars:      cpro_vscalt kg/(m s)

  \warning: if the scalar is the temperature, cpro_vscalt corresponds
  to its conductivity (Lambda) in W/(m K)

  The types of boundary faces at the previous time step are available
    (except at the first time step, where arrays \c itypfb and \c itrifb have
    not been initialized yet)

  It is recommended to keep only the minimum necessary in this file
    (i.e. remove all unused example code)


  \section example1_comp Molecular viscosity varying with temperature

  The values of the molecular viscosity are provided as a function of
  the temperature. All variables are evaluated at the cell centers.

  Here is the corresponding code:

  \snippet cs_user_physical_properties-compressible_flow.f90 example_1


  \section example2_comp Molecular volumetric viscosity varying with temperature

  The values of the molecular volumetric viscosity are provided as a function
  of the temperature. All variables are evaluated at the cell centers.

  Here is the corresponding code:

  \snippet cs_user_physical_properties-compressible_flow.f90 example_2


  \section example3_comp Isobaric specific heat varying with temperature

  The values of the isobaric specific heat values are provided as a function
  of the temperature. All variables are evaluated at the cell centers.

  \warning:
  do not discard the call to the subroutine 'usthht' at the end of this
  example: its purpose is to calculate the isochoric specific heat.
  Indeed, this variable needs to be computed from the isobaric specific heat
  using the thermodynamics laws.

  Here is the corresponding code:

  \snippet cs_user_physical_properties-compressible_flow.f90 example_3


  \section example4_comp Molecular thermal conductivity varying with temperature

  The values of the molecular thermal conductivity are provided as a function
  of the temperature. All variables are evaluated at the cell centers.

  Here is the corresponding code:

  \snippet cs_user_physical_properties-compressible_flow.f90 example_4


  \section example5_comp Molecular diffusivity of user-defined scalars varying with temperature

  The molecular diffusivity can be set for all the user-defined scalars
  <b>except</b>:
    - temperature and enthalpy (already dealt with above: for these
      variables, the 'diffusivity' is the thermal conductivity)
    - variances of the fluctuations of another scalar variable (the
      diffusivity is assumed to be equal to that of the associated
      scalar)
  The values of the molecular diffusivity are provided as a function
  of the temperature. All variables are evaluated at the cell centers.

  Here is the corresponding code:

  \snippet cs_user_physical_properties-compressible_flow.f90 example_5


*/
