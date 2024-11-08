<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
-->

\page cs_ug_cdo_solidification Solidification modelling or melting of solid with CDO schemes

[TOC]

<!--
    References used in this page
-->

[VolPra]: https://www.sciencedirect.com/science/article/abs/pii/0017931087903176
[BF23]: https://hal.science/hal-04028238

Introduction
============

Modelling phase changes between a solid and liquid phase is possible
with code_saturne. Historically, this module has been developped for
modelling the solidification process, hence the name "solidification"
module. Nevertheless, the melting of a solid is also possible.

The solidification module described in this page relies on CDO schemes
and especially, CDO face-based schemes.

To set-up a computation taking into account the solidification
process, one has to update the \ref cs_user_parameters.c user source
file and edit at least the user-defined functions \ref cs_user_model
and \ref cs_user_finalize_setup in simple cases. In more complex
cases, editing the function \ref cs_user_parameters and defining
functions to describe the boundary conditions or the material
properties should be necessary.

Model of solidification and its options {#cs_ug_cdo_solidification_model}
=======================================

Several models are available according to the type of phenomena at stake.
Here are listed the available solidification models (first parameter of the
function \ref cs_solidification_activate):

- \ref CS_SOLIDIFICATION_MODEL_STEFAN : a pure thermal model,
  i.e. without a fluid motion. A rough phase transition is considered
  (the solidus temperature is equal to the liquidus temperature).

- \ref CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 : a model involving
  momentum, mass and energy conservation equations. A Boussinesq
  approximation is used for the Navier--Stokes model. The thermal
  source term is linearized. A pure component is considered (not seen
  as an alloy). The solidification process hinges on a head loss
  described as a Darcy-like source term in the momentum equation (see
  the [Voller and Prakash article](https://www.sciencedirect.com/science/article/abs/pii/0017931087903176)
  for more details).

- \ref CS_SOLIDIFICATION_MODEL_VOLLER_NL : This is a variation of the
  \ref CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 model with a
  non-linear thermal source term (a Picard algorithm is used to update
  the source term).

- \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY : This model hinges on the
  Voller and Prakash model for the dynamics but adds an additional
  equation for the solute transport since the mixture is assumed to be
  a binary alloy. The solute concentration has also an effect on the
  Boussinesq term. The phase diagram is more complex since a binary
  alloy is considered. One also handles the eutectic transformation
  when the solidification path reaches the eutectic plateau.

Model options {#cs_ug_cdo_solidification_model_options}
-------------

Besides the main solidification model, one can specify several
modelling options. Be aware that some options are only relevant when a
given model is set. Here are listed the available options (second
parameter of the function \ref cs_solidification_activate):

- \ref CS_SOLIDIFICATION_NO_VELOCITY_FIELD : This flag disables the
  resolution of the Navier-Stokes equations. If used with the
  \ref CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 and
  \ref CS_SOLIDIFICATION_MODEL_VOLLER_NL models, this corresponds to a
  Stefan-like model with a phase transition having a different value
  for the solidus and liquidus temperatures.

- \ref CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM (only make sense when the
  model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY is activated). The bulk
  concentration is treated with a source term related to an explicit
  advection of the quantity \f$ (C_m - C_l) \f$.

Automatic post-processing {#cs_ug_cdo_solidification_post}
-------------------------

Here are listed the available quantities which can be automatically
post-processed (third parameter of the function \ref cs_solidification_activate):

- \ref CS_SOLIDIFICATION_POST_CELL_STATE
- \ref CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE
- \ref CS_SOLIDIFICATION_POST_ENTHALPY

The following options make sense only when the
\ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY is activated

- \ref CS_SOLIDIFICATION_POST_CBULK_ADIM
- \ref CS_SOLIDIFICATION_POST_CLIQ
- \ref CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE
- \ref CS_SOLIDIFICATION_POST_SEGREGATION_INDEX


Activate the solidification module {#cs_ug_cdo_solidification_activate}
==================================

The first step when settings a solidification model relying on CDO
schemes is to set the CDO mode to \ref CS_PARAM_CDO_MODE_ONLY

\snippet cs_user_parameters-cdo-condif.cpp param_cdo_activation

This is done at the begining of the function \ref cs_user_model.

Then, the second step is to add a call to the function
\ref cs_solidification_activate in order to activate the CDO
solidification module. There are eight parameters to set :

1. The type of solidification model to consider. This implies which
   equations are to be solved and which variables are
   involved. Choices are detailed in \ref cs_solidification_model_t

2. An optional flag to specify advanced user-defined options (set 0 if
   nothing has to be added);

3. An optional flag to specify automatic post-processing dedicated to
   this module (set 0 if nothing has to be added);

4. The domain boundaries (see \ref cs_boundary_t for more details);

5. The Navier-Stokes model (see \ref cs_navsto_param_model_t for more
   details). The default choice should be
   \ref CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES but a Stokes model
   (set \ref CS_NAVSTO_MODEL_STOKES) is also possible if the advection
   term can be neglected.

6. Additional user-defined options for the Navier-Stokes model (by
   default the Boussinesq approximation is set even if 0 is given as
   parameter). The recommended value for this parameter is 0.

7. The algorithm used to couple velocity and pressure
   fields. Available choices are detailed in \ref cs_navsto_param_coupling_t

8. An optional flag to specify automatic post-processing for the
   Navier-Stokes module. Available options are detailed in
   \ref cs_navsto_param_post_bit_t).

Here are several examples of activation of the CDO solidification
module according to the model and its related options.

**Ex. 1** The first example describes how to set a "binary alloy"
model with several related options.

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_activate_solidification_binary

**Ex. 2** The second example describes how to set a "Voller" model
with several related options.

  \snippet cs_user_parameters-cdo-solidification.cpp param_cdo_activate_solidification_voller

**Ex. 3** The third example describes how to set a "Voller" model
without a velocity field in the "liquid" area.

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_activate_solidification_voller_no_vel


Model settings {#cs_ug_cdo_solidification_set_model}
==============

Set a strategy {#cs_ug_cdo_solidification_strategy}
--------------

Excepted for the Stefan model, it is possible to change the strategy
to update quantities such as the way to compute the new thermal source
term.  There are three strategies available (\ref cs_solidification_strategy_t
for more details). The default strategy is \ref CS_SOLIDIFICATION_STRATEGY_LEGACY

- \ref CS_SOLIDIFICATION_STRATEGY_LEGACY : Original strategy considered in
  the Finite Volume scheme

- \ref CS_SOLIDIFICATION_STRATEGY_TAYLOR : Update of the thermal source term
  or the liquid fraction is performed using a Taylor expansion. Only
  available with a \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY

- \ref CS_SOLIDIFICATION_STRATEGY_PATH : Update of the thermal source (and of
  the liquid fraction in case of a \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY)
  relying on the knowledge of the current and previous state.

Here is an example to set the strategy:

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_set_strategy


Stefan model {#cs_ug_cdo_solidification_set_stefan}
------------

  The Stefan model handles the solidification process with the approximation
  that the liquidus and solidus temperature are the same. Thus, the transition
  from the solid to liquid or from the liquid to solid is rough.

Voller and Prakash model {#cs_ug_cdo_solidification_set_voller}
------------------------

The main function to set a Voller and Prakash model is \ref
cs_solidification_set_voller_model (follow the link to get a
description of the parameters). This function is mandatory and has to be
located in \ref cs_user_model after calling \ref cs_solidification_activate

- Parameters 1 and 2 correspond to the settings of the Boussinesq term (two
  coefficients related to a reference temperature and a thermal dilation
  coefficient in \f$ K^-1\ \f$).

- Parameters 3 and 4 are the solidus and liquidus temperature describing the
  behavior of the liquid fraction with respect to the temperature.

- Parameters 5 and 6 correspond to the latent heat of the alloy (assumed to
  be constant) and the secondary dendrite arm spacing (a model parameter
  taking into account micro-segregation phenomena).

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_set_voller

Please notice that a simplified version of the function exists in case of
purely thermal model (activated when the flag \ref
CS_SOLIDIFICATION_NO_VELOCITY_FIELD has been set).

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_set_voller_no_vel

Non-linear Voller and Prakash model {#cs_ug_cdo_solidification_set_voller_nl}
-----------------------------------

The setting of the non-linear variant of the Voller and Prakash model relies
on the main function as the linear one (i.e. \ref
cs_solidification_set_voller_model or \ref
cs_solidification_set_voller_model_no_velocity in a more specific situation).

### Advanced usage

To access more settings, it is possible to retrieve the structure managing
the voller model and to specify advanced parameter settings in order to set
the non-linear iterative algorithm for instance.

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_nl_voller_advanced


Binary alloy model {#cs_ug_cdo_solidification_set_binary}
------------------

The main function to set a binary alloy model is
\ref cs_solidification_set_binary_alloy_model (follow the link to get a
description of the parameters). This function is mandatory and has to be
located inside \ref cs_user_model after calling \ref cs_solidification_activate

- Parameters 1 and 2 correspond to the name of the equation related to the
  solute transport and the name of the unknow associated to this equation
  (this will be the name of the variable field).

- Parameters 3 to 6 correspond to the settings of the Boussinesq term (two
  contributions: the classical thermal one given by a reference temperature
  and a thermal dilation coefficient in \f$ K^-1 \f$; the solutal
  contribution to the Boussinesq term given by a reference concentration and
  a solutal dilation coefficient)

- Parameters 7 to 10 correspond to the main parameters describing the phase
  diagram (the partition coefficient between solid and liquid phase, the
  liquidus slope, the eutectic temperature and the melting temperature for a
  pure material which corresponds to the one obtained when one of the
  component of the binary alloy is not present).

- Parameter 11 is the value of the diffusivity for the solute transport
  equation

- Parameters 12 and 13 correspond to the latent heat of the alloy (assumed
  to be constant) and the secondary dendrite arm spacing (a model parameter
  taking into account micro-segregation phenomena).

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_set_binary_alloy

###  Advanced usage {#cs_ug_cdo_solidification_set_binary_x}

To access more settings, it is possible to retrieve the structure
managing the binary alloy model and to specify advanced parameter
settings.

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_binary_advanced



Property settings {#cs_ug_cdo_solidification_set_property}
=================

The creation of all properties used the solidification model is made
automatically. The definition is then performed in the function \ref
cs_user_finalize_setup

In simple cases, the main properties are constant (uniform and steady-state)
and the definition is thus simply done in the example below. Please notice
that predefined properties have a macro to store their name. Here are used:

 \ref CS_PROPERTY_MASS_DENSITY for the (volumetric) mass density which
 corresponds to the property named "mass_density" (this is different from the
 legacy Finite Volume where "density" is used);

 \ref CS_NAVSTO_LAM_VISCOSITY for the laminar (dynamic) viscosity which
 corresponds to the property named "laminar_viscosity";

 \ref CS_THERMAL_CP_NAME for the thermal capacity which corresponds to the
 property named "thermal_capacity";

 \ref CS_THERMAL_LAMBDA_NAME for the thermal conductivity which corresponds
 to the property named "thermal_conductivity".

snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_properties


Equation settings {#cs_ug_cdo_solidification_set_eq}
=================

The last step is to set the initial condition and the boundary
conditions. By default, all variable fields are initially set to
zero. For the boundary conditions, if nothing is set, then the default
boundary condition is applied. For the thermal equation, the default
boundary condition is a no flux (i.e. a homogeneous Neumann boundary
condition).

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_thermal_eq

When a solute transport equation is added (this done automatically when
calling the function \ref cs_solidification_set_binary_alloy_model), the
default boundary condition for this equation is a no flux condition as well.

Here is an example how to set the initial solute concentration in the domain.

\snippet cs_user_parameters-cdo-solidification.cpp param_cdo_solidification_solute_eq

For the Navier-Stokes equation, the default boundary condition is defined when
the domain boundaries are defined in the function \ref cs_user_model using the
function \ref cs_boundary_set_default There are two possibilities:

- \ref CS_BOUNDARY_SYMMETRY
- \ref CS_BOUNDARY_WALL


To go further
=============

A more detailed description of the model and a comparison between the
Finite Volume and CDO approaches is available [here](https://hal.science/hal-04028238)
