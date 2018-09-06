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
  \page cavit Data setting for the cavitation module

  \section cavitation_h_intro Introduction

  The cavitation module is based on a homogeneous mixture model.
  The physical properties (density and dynamic viscosity) of the mixture depend
  on a resolved void fraction and constant reference properties of the liquid
  phase and the gas phase.
  The void fraction is given by an advection equation with a
  vaporization/condensation source/sink term. This term is modeled by the
  Merkle's model. The module also integrates the eddy-viscosity correction of
  Reboud.

  \section cavit_activ Activation of the module

  The module can be activated in the \ref usipph routine in
  \ref cs_user_parameters.f90. The corresponding keyword is icavit in the
  \ref optcal module.
  This keyword can take the values:
   - \ref optcal::icavit "icavit" = -1: module desactivated (single-phase flow).
   - \ref optcal::icavit "icavit" =  0: the module is activated but there is no
        vaporization/condensation source term. The void fraction is only
        advected by the mixture velocity.
   - \ref optcal::icavit "icavit" =  1: the module is activated and the Merkle
        vaporisation/condensation source/sink term is taken into account.

  \section cavit_parameters Cavitation module specific parameters.

  When the module is activated, its specific input parameters should be set in
  the \ref usipsu routine of \ref cs_user_parameters.f90 file. An example is given
  in cs_user_parameters-cavitation.f90.

  \subsection cavit_phprop Homogeneous mixture physical properties

  As soon as icavit  \f$ \ge 0 \f$, the reference density, in \f$ kg/m^3\f$ , and molecular viscosity, \f$ kg/(m\cdot s)\f$, of the liquid phase and the gas phase should be set. For instance:

  \snippet cs_user_parameters-cavitation.f90 phprop_l
for the liquid and:
  \snippet cs_user_parameters-cavitation.f90 phprop_g
for the gas phase.

  \subsection cavit_source Model parameters of the vaporization term (Merkle model)

  When icavit = 1, Merkle's model parameters should be set.
  The Merkle model is base on a barotropic law for the density (see \ref cavitation.f90). In that way, its principal parameter is the saturation pressure of the fluid, in \f$ kg/(m\cdot s^2)\f$. For instance, the saturation pressure of the water at twenty celcius degrees is:

  \snippet cs_user_parameters-cavitation.f90 presat

  Merkle's model also requires a reference length scale and velocity of the flow. For instance:

  \snippet cs_user_parameters-cavitation.f90 scales_inf

  These scales are integral scales. For instance, considering the cavitating flow across a foil in a duct, the reference velocity should be the bulk velocity and the reference length scale should be the chord of the foil.

  \subsection cavit_turb Interaction with turbulence

  As soon as icavit \f$ \ge 0 \f$, the mixture eddy-viscosity correction proposed by Reboud can be accounted for:

  \snippet cs_user_parameters-cavitation.f90 reboud_activ

If icvevm = 0, the Reboud correction is deactivated. Using an eddy-viscosity model (see \ref turbulence), this option is recommended, such that icvevm = 1 is the default setting. Of course, this option has no effect for second moment closure or large eddy simulations. Note that the the coefficent mcav of the reboud correction (see \ref cavitation.f90) can also be adjust in the \ref usipsu routine.

  \subsection cavit_numerics Numerical options

  Advanced numerical parameters may also be set in this routine, if necessary. The concerned variables are listed in \ref cav_numerics.

*/
