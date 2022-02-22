<!--
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
-->

\page cs_gwf_cdo Groundwater flow module using CDO schemes

[TOC]

<!--
    References used in this page
-->

[Bonel14]: https://hal.archives-ouvertes.fr/tel-01116527
[BoFoM18]: https://dx.doi.org/10.1016/j.compfluid.2018.03.026
[BoFoM18_hal]: https://hal.archives-ouvertes.fr/hal-01624854

Introduction {#sec_cdo_gwf_intro}
============

The Hydrogeology module of code_saturne is a numerical model for groundwater
flow and solute transport in continuous porous media. The flow part is based on
the Richards equation, derived from the Darcy law and the conservation of
mass. The transport part is based on the the classical advection-diffusion
equation of tracers, slightly modified to account the specificities of
groundwater transport.

This module can be used to simulate transfers of water and solutes in several
saturated and/or unsaturated porous media. The flow part can be steady or
unsteady, with isotropic or anisotropic permeabilities and allows for any type
of soil water retention model thanks to a user-defined model. Two classical
models are predefined: the saturated model and the van Genuchten-Mualen
model. The transport part considers dispersion, sorption and radioactive
decay. The partition between soil and water phases can be modeled by a
classical Kd approach model. Additionaly solute precipitation/dissolution
phenomena can also be taken into account by an instantaneous model.

Physical concepts and equations are presented in the [theory guide](../../theory.pdf)

The groundwater flow module (GWF) relies on CDO vertex-based or CDO
vertex+cell-based discretization schemes. Here is listed a set of references
useful for the reader willing to get a better understanding of the mathematical
concepts underpinning CDO schemes.

* [**PhD thesis**: New Polyhedral Discretisation Methods applied to the Richards Equation: CDO Schemes in Code Saturne][Bonel14]

* [**Article:** New Polyhedral Discretisation Methods applied to the Richards Equation: CDO Schemes in Code Saturne][BoFoM18] ([**HAL** preprint version][BoFoM18_hal])


To set-up a GWF computation, one has to update the cs_user_parameters.c file
and edit the function \ref cs_user_model at least in simple cases. In more
complex cases, editing \ref cs_user_finalize_setup should be necessary.


Activate the GWF module {#sec_cdo_gwf_activate}
=======================

The first step is to activate the CDO module in the function \ref cs_user_model
(please refer to \ref cs_user_parameters_h_cdo_activation).

Then, one has to activate the groundwater flow (GWF) module in the function
\ref cs_user_model.  The function to call is \ref cs_gwf_activate.

There are three parameters:

1. The main hydraulic model to consider (i.e. which equations have to be
   solved). Please refer to \ref cs_gwf_model_type_t. There are currently two
   models available :

  * \ref CS_GWF_MODEL_SATURATED_SINGLE_PHASE : The simplest model (the Richards
    equation becomes a simple steady diffusion equation). The name of the
    equation which is automatically created is "Richards".

  * \ref CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE . The name of the equation which
    is automatically created is "Richards".

2. An optional flag to specify a physical phenomena to take into account or to
   specify a numerical treatment to apply. If no option is needed, then simply
   set 0. Here are listed the available option flags:

  - \ref CS_GWF_GRAVITATION
  - \ref CS_GWF_FORCE_RICHARDS_ITERATIONS
  - \ref CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE
  - \ref CS_GWF_ENFORCE_DIVERGENCE_FREE

3. An optional flag to specify the activation of automatic postprocessing for
   quantities of interest related to groundwater flows. Set 0 if no additional
   postprocessing is requested. Here are listed the available options:

  - \ref CS_GWF_POST_SOIL_CAPACITY
  - \ref CS_GWF_POST_LIQUID_SATURATION
  - \ref CS_GWF_POST_PERMEABILITY
  - \ref CS_GWF_POST_DARCY_FLUX_BALANCE
  - \ref CS_GWF_POST_DARCY_FLUX_DIVERGENCE
  - \ref CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY

___________

**Remark:** If a \ref CS_GWF_MODEL_SATURATED_SINGLE_PHASE is set at the
  activation step, then one expects that all soil models are defined by the
  type \ref CS_GWF_SOIL_SATURATED
___________


Examples of activation of the GWF module
----------------------------------------

**Example 1:** _Activate the GWF model with a fully saturated single-phase
flow model and no other option._

\snippet cs_user_parameters-cdo-gwf.c param_cdo_activate_gwf

**Example 2:** _Second example: Activate the GWF model with an unsaturated
single-phase flow model without any additional option._

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_activate_gwf_uspf

**Example 3:** _Activate the GWF model with an unsaturated single-phase flow
model. Moreover, one takes into account the gravity effect and the
postprocessing of the soil permeability as well as the Darcy flux across domain
boundaries._

\snippet cs_user_parameters-cdo-gwf.c param_cdo_activate_gwf_b




Soils
=====

The second step is to **add at least one new soil**. The add of soil(s) should
be done before adding tracers inside the function \ref cs_user_model and after
having activated the GWF module. Two functions are available to add a new soil:

* [cs_gwf_add_iso_soil](@ref cs_gwf_add_iso_soil) when the soil is modelled by
  an isotropic absolute (or intrinsic) permeability

* \ref cs_gwf_add_aniso_soil when the soil is modelled by an anisotropic
  absolute (or intrinsic) permeability (_i.e._ one has to specific a 3x3
  tensor)

These two functions have a similar set of parameters:

1. The name of the volume zone associated to this soil (to add a volume, one
   can either use the GUI or use \ref cs_volume_zone_define inside \ref
   cs_user_zones ; see \ref cs_user_zones_volume_simple for more details)

2. The value of the bulk mass density of the soil. This is only useful when a
   tracer is considered. If there is no tracer, one can set `1` for instance.

3. The value of the absolute permeability. In case of an isotropic soil, this
   is a scalar and in case of an anisotropic soil, one expects a tensor.

4. The value of the soil porosity (which is equivalent to the saturated
   moisture or the max. liquid saturation in single-phase flows)

5. The type of soil model. There are two predefined soil models and one
   user-defined soil model

  - \ref CS_GWF_SOIL_SATURATED
  - \ref CS_GWF_SOIL_GENUCHTEN
  - \ref CS_GWF_SOIL_USER


Examples of settings for a predefined soil model
------------------------------------------------

### Case of a saturated model

The saturated model is the simplest model.

**Example 1:** _A saturated soils defined by an isotropic permeability on all
the computational domain._

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_add_define_iso_saturated_soil

**Example 2:** _Two saturated soils defined by an anisotropic (saturated)
permeability_

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_add_define_aniso_saturated_soil


### Case of a Van Genuchten-Mualen model

Soils which behave according to a Van Genuchten-Mualen model are specified in
two steps: a call to \ref cs_gwf_add_iso_soil and then a call to \ref
cs_gwf_soil_set_genuchten_param to specifiy the parameters associated to this
model.

**Example 3:** _Soil relying on a Van Genuchten-Mualen and considering a
isotropic permeability_

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_add_define_genuchten_soil


Example of settings for a user-defined soil
-------------------------------------------

If the predefined models are not sufficient, it is possible to add a
user-defined soil. In this case, the add of the new soil is made as follows

1. Define a structure to get an access to the parameters defining the soil model

2. Add a new soil (\ref cs_gwf_add_iso_soil or \ref cs_gwf_add_aniso_soil )

3. Call \ref cs_gwf_soil_set_user to specify the structure to use, the function
   to update quantities related to the hydraulic model and if needed a function
   to free the structure storing the soil parameters.

These three steps are performed inside the function \ref cs_user_model

Here is a complete example of a user-defined soil model (called hereafter Tracy
since it has been designed by F. T. Tracy in [this
article](https://doi.org/10.1016/0022-1694(94)02674-Z)).

### Define a structure to store the model parameters

_Example of the structure used to handle the soil model parameters_

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_gwf_tracy_struct

### Add a user-defined soil

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_gwf_add_user_soil

with the two requested functions (defined for instance as a static function in
the file cs_user_parameters.c). These functions have to fullfill the prototype
defined in \ref cs_gwf_soil_update_t (for the update of the soil properties)
and in \ref cs_gwf_soil_free_param_t (for the free of the soil parameter
structure).

Here is an example of how to update soil properties (function called
_tracy_update_)

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_gwf_set_user_update_soil

and an example of how to free the soil parameter structure (function called
_tracy_free_param_)

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_gwf_set_user_free_soil

### Further settings (initial and boundary conditions)

In this example, we also give some highlights on how this soil structure can be
used to further set the problem for instance to specify the initial and
boundary conditions.
This step is made in the function \ref cs_user_finalize_setup

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_gwf_get_user_soil

where the two functions used to define either the boundary condition (called
"get_bc") in the example or the initial condition (called "get_ic") in the
example follow a predefined prototype (see \ref cs_analytic_func_t)

Here are collected two examples of such functions:

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_gwf_set_bc_analytic

\snippet cs_user_parameters-cdo-gwf_user_soil.c param_cdo_set_ic_by_analytic



Tacers
======

The third step (which is not mandatory) is to add tracer(s) thanks to the
function \ref cs_gwf_add_tracer This tracer will be advected by the Darcy flux
arising from the Richards equation.

There are currently two models :

  - a default model (the predefined one; see \ref cs_gwf_cdo_predef_tracer)
  - a user-defined model (see \ref cs_gwf_cdo_user_tracer)

The first parameter in \ref cs_gwf_add_tracer is a flag which can be built with

  - \ref CS_GWF_TRACER_USER (to switch to a user-defined tracer)
  - \ref CS_GWF_TRACER_PRECIPITATION (to add the precipitation effect)


Predefined tracers
------------------

Here is a simple example for a standard tracer which can be added in the
function \ref cs_user_model

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_add_tracer

Remark: Get a tracer structure.

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_get_tracer


User-defined tracers
--------------------

TODO




Automatic postprocessings
=========================

It is possible to activate or add an automatic post-processing of several
quantities of interest related to groundwater flows. Here are available flags
to activate through the usage of \ref cs_gwf_set_post_options

  - \ref CS_GWF_POST_SOIL_CAPACITY
  - \ref CS_GWF_POST_LIQUID_SATURATION
  - \ref CS_GWF_POST_PERMEABILITY
  - \ref CS_GWF_POST_DARCY_FLUX_BALANCE
  - \ref CS_GWF_POST_DARCY_FLUX_DIVERGENCE
  - \ref CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY

\snippet cs_user_parameters-cdo-gwf.c param_cdo_post_gwf



Helper functions {#sec_cdo_gwf_helper}
================

Helper functions for soils {#sec_cdo_gwf_helper_soil}
--------------------------

Get a soil structure from its name.

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_get_soil

There is a similar which retrieve the soil structure from its id (see \ref
cs_gwf_soil_by_id).

Helper functions for tracers {#sec_cdo_gwf_helper_tracer}
----------------------------

Get a tracer structure from its name.

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_get_tracer
