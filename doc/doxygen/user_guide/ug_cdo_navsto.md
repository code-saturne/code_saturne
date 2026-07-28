<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

\page cs_ug_cdo_navsto Using CDO schemes for solving the Navier-Stokes system

[TOC]

<!--
    References used in this page
-->

[Mila20]: https://hal.science/tel-03080530v2
[MiBoE22]: https://hal.science/hal-03215118v1

Introduction {#sec_ug_cdo_navsto_intro}
============

The Navier-Stokes equations can be solved using **CDO face-based**
discretizations. Other CDO/HHO spatial discretizations are not available.
Currently, there is no turbulence modeling natively supported when using CDO
schemes. This is a work in progress; however, turbulence modeling can be
achieved by coupling CDO face-based schemes with legacy Finite Volume (FV)
schemes for the turbulence equations.

The rationale for setting up Navier-Stokes computations is very similar to the
one explained for user-defined equations [here](@ref cs_ug_cdo_hho_base). This
guide focuses solely on the specific features arising from the Navier-Stokes
equations.


Settings done in cs_user_model()
===================

The NavSto module is activated using the function \ref
cs_navsto_system_activate. The choice of parameters in this function is crucial.
For instance, here are two examples demonstrating how to activate and configure
the main the physical and numerical parameters of the NavSto module:

\snippet cs_user_parameters-cdo-navsto.cpp param_cdo_navsto_activate

or

\snippet cs_user_parameters-cdo-navsto.cpp param_cdo_navsto_activate2

The first parameter is the structure managing the domain boundaries. An example
of the domain boundary settings is described [here](@ref ug_cdo_sec_domain_boundaries).


### Choosing your model and coupling algorithm {#sec_ug_cdo_navsto_model_choice}

The two most important parameters are the physical model (\ref cs_navsto_param_model_t) and the
velocity-pressure coupling algorithm (\ref cs_navsto_param_coupling_t).

The second parameter of the function `cs_navsto_system_activate` specifies the physical model.

| Model (`cs_navsto_param_model_t`) | Description | When to use? |
| :--- | :--- | :--- |
| `CS_NAVSTO_MODEL_STOKES`                      | Stokes flow (Re << 1). Inertia terms are neglected. | Very viscous flows, microfluidics. Linear problem. |
| `CS_NAVSTO_MODEL_OSEEN`                       | Linearization of Navier-Stokes equations around a given velocity field. | Used in fixed-point algorithms. |
| `CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES` | Full incompressible Navier-Stokes equations. Non-linear. | General purpose for incompressible flows. |

The *base* model can be updated using the third parameter, which is a flag
built from the following bitwise options (\ref cs_navsto_param_model_bit_t):
meaning that these options can be combined:
    - \ref CS_NAVSTO_MODEL_STEADY to specify that the model is steady-state (by
      default, it is unsteady).
    - \ref CS_NAVSTO_MODEL_GRAVITY_EFFECTS
    - \ref CS_NAVSTO_MODEL_BOUSSINESQ

Other advanced settings are detailed in \ref cs_navsto_param_model_bit_t.

The fourth parameter of the function `cs_navsto_system_activate` specifies the
velocity-pressure coupling algorithm (\ref cs_navsto_param_coupling_t). The
options are:

| Coupling (`cs_navsto_param_coupling_t`) | Principle | Pros / Cons | Comments |
| :--- | :--- | :--- | :--- |
| \ref CS_NAVSTO_COUPLING_MONOLITHIC | Velocity and pressure are solved simultaneously in a single large linear system (saddle-point).    | **+** Robust, strong coupling. <br> **-** The linear system is complex and memory-intensive. Requires specific preconditioners. | More details and a mathematical analysis are gathered in \cite MiBoE22 and \cite Mila20 |
| \ref CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY | An artificial time derivative is added to the continuity equation to couple velocity and pressure. | **+** Transforms the problem into a simpler transport-diffusion-like system. <br> **-** Mainly for steady-state cases. The artificial compressibility parameter must be chosen carefully. | More details and a mathematical analysis are gathered in \cite MiBoE22 and \cite Mila20 |
| \ref CS_NAVSTO_COUPLING_PROJECTION | | | **Deprecated**. Do not use anymore. `AFS` solver for saddle-point problem is a better choice |



The last parameter specifies predefined post-processing operations. Similar to
the third parameter, this is a flag constructed from the following bitwise
options (\ref cs_navsto_param_post_bit_t):
    - \ref CS_NAVSTO_POST_VELOCITY_DIVERGENCE
    - \ref CS_NAVSTO_POST_KINETIC_ENERGY
    - \ref CS_NAVSTO_POST_VORTICITY
    - \ref CS_NAVSTO_POST_VELOCITY_GRADIENT
    - \ref CS_NAVSTO_POST_STREAM_FUNCTION (This adds an equation named
      *streamfunction_eq* and its variable field named *stream_function*. This
      equation relies on a CDO vertex-based discretization.)
    - \ref CS_NAVSTO_POST_HELICITY
    - \ref CS_NAVSTO_POST_ENSTROPHY
    - \ref CS_NAVSTO_POST_MASS_DENSITY
    - \ref CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE
    - \ref CS_NAVSTO_POST_PRESSURE_GRADIENT
    - \ref CS_NAVSTO_POST_BOUNDARY_STRESS


Further settings for the numerics are also available using \ref
cs_navsto_param_num_bit_t.

Predefined equations associated with the Navier-Stokes equations are:
    - The *momentum* equation, which is automatically added when the module is
    activated with either a **monolithic** or **artificial compressibility**
    velocity-pressure coupling.



In all cases, a vector-valued field named _velocity_ and a scalar-valued field
named _"pressure"_ are created. Moreover, several properties are added:
    - the property _mass_density_ (a \ref cs_property_t structure);
    - the property _laminar viscosity_ (a \ref cs_property_t structure);
    - the properties _turbulent_viscosity_ and _total_viscosity_ are added
      if the turbulence model is different from the laminar one (cf. \ref cs_turb_model_t);
    - the property _graddiv_coef_ (a \ref cs_property_t structure) when
      **artificial compressibility** is set;

along with the advection field _mass_flux_ (a \ref cs_adv_field_t structure).


Settings done in cs_user_finalize_setup()
===================

Boundary conditions
-------------------

For the NavSto module, boundary conditions are set as follows (direct access
to the equation is not required, as it depends on the chosen velocity-pressure
coupling algorithm):

\snippet cs_user_parameters-cdo-navsto.cpp param_cdo_navsto_bc1

where the function `_vel_def` is defined as follows:

\snippet cs_user_parameters-cdo-navsto.cpp param_cdo_navsto_vel_function

Definition of source terms
---------------------

Since the equation to which this source term is applied depends on the choice
of the velocity-pressure coupling algorithm, the procedure differs slightly
from the one used for user-defined equations:

\snippet cs_user_parameters-cdo-navsto.cpp param_cdo_navsto_st1

where the function `_src_def` is defined as follows:

\snippet cs_user_parameters-cdo-navsto.cpp param_cdo_navsto_st_function


Settings done in cs_user_parameters() {#cs_ug_cdo_navsto_set_param}
===================

The setup procedure is similar to the one detailed in [this section](@ref cs_ug_cdo_hho_base_user_param).

This section lists some aspects specific to the NavSto module.

Set the strategy to solve the Navier-Stokes system when a monolithic coupling is used
-------------------

When a *monolithic* velocity-pressure coupling is chosen, the resulting linear
system to solve is a *saddle-point* problem. This class of linear systems
requires specific choices of preconditioners and solvers. The default settings
are not always optimal for performance. The configuration of saddle-point
problems is detailed [here](@ref cs_ug_cdo_sles_saddle).


To go further
=============

The detailed description of CDO schemes, their mathematical analysis, and
numerical benchmark results are available in the following PhD thesis:

* [**PhD thesis**: "Compatible Discrete Operator schemes for the unsteady incompressible Navier-Stokes equations"][Mila20] \cite Mila20


Additional publications:
* [**Article**: "Artificial compressibility methods for the incompressible Navier-Stokes equations using lowest-order face-based schemes on polytopal meshes"][MiBoE22] \cite MiBoE22
