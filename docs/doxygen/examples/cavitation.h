/*============================================================================
 * code_saturne documentation page
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
  \page cavit Data setting for the cavitation model

  \section cavitation_h_intro Introduction

  The cavitation model is based on a homogeneous mixture model, and a sub-model
  of the Volume of Fluid model.
  The physical properties (density and dynamic viscosity) of the mixture depend
  on a resolved void fraction and constant reference properties of the liquid
  phase and the gas phase.
  The void fraction is given by an advection equation with a
  vaporization/condensation source/sink term. This term is modeled by the
  Merkle's model. The model also integrates the eddy-viscosity correction of
  Reboud.

  \section cavit_activ Activation of the model

  The module can be activated in \ref cs_user_model in
  \ref cs_user_parameters.c as show below:
  \snippet cs_user_parameters-base.c enable_cavit

  Data structure is defined in \ref vof, \ref vof_mixture_properties and
  \ref cavitation.
  Cavitation modelling main feature consists in modelling
  vaporisation/condensation with Merkle model providing source / sink term for
  the void fraction equation.

  \section cavit_parameters Cavitation module specific parameters.

  When cavitation model is enabled, specific input parameters can be set in
  \ref cs_user_parameters in \ref cs_user_parameters.c file as shown below:

  \subsection cavit_phprop Homogeneous mixture physical properties

  The reference density, in \f$ kg/m^3\f$ , and molecular viscosity, \f$ kg/(m\cdot s)\f$,
  of the liquid phase and the gas phase should be set. For instance:

  \snippet cs_user_parameters-base.c phprop

  Other parameters, specific to cavitation are stored in \ref cs_cavitation_parameters_t
  structure. A pointer to this structure should be retrieved as follows:

  \snippet cs_user_parameters-base.c cavit_param

  \subsection cavit_source Model parameters of the vaporization term (Merkle model)

  Merkle's model parameters should be set.
  Merkle's model is based on a barotropic law for the density (see \ref cavitation.f90).
  In that way, its principal parameter is the saturation pressure of the fluid,
  in \f$ kg/(m\cdot s^2)\f$. For instance, the saturation pressure of the water at twenty
  celcius degrees is set below:

  \snippet cs_user_parameters-base.c presat

  Merkle's model also requires a reference length scale and velocity of the flow.
  For instance:

  \snippet cs_user_parameters-base.c scales_inf

  These scales are integral scales. For instance, considering the cavitating
  flow across a foil in a duct, the reference velocity should be the bulk
  velocity and the reference length scale should be the chord of the foil.

  \subsection cavit_turb Interaction with turbulence

  The mixture eddy-viscosity correction proposed by Reboud can be accounted for
  as shown below:

  \snippet cs_user_parameters-base.c reboud_activ

  Using an eddy-viscosity model (see \ref turbulence), this option is
  recommended and is hence a default setting. Of course, this option has no
  effect for second moment closure or large eddy simulations. Note that the
  coefficent mcav of the reboud correction (see \ref cavitation.f90) can also
  be adjusted in the \ref cs_user_parameters function.

  \subsection cavit_numerics Numerical options

  Advanced numerical parameters may also be set in this function, if necessary.
  The concerned variables are listed in \ref cav_numerics.

*/
