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

\page cs_ug_cdo_navsto Using CDO schemes for solving the Navier-Stokes system

[TOC]

<!--
    References used in this page
-->

[Mila20]: https://hal.science/tel-03080530v2
[MiBoE22]: https://hal.science/hal-03215118v1

Introduction {#sec_ug_cdo_navsto_intro}
============

The Navier--Stokes equations can be solved using **CDO face-based**
discretizations. Other space discretizations are not available. Up to
now, there is **no turbulence modelling** when using CDO schemes. This
is a work in progress.

The rationale to set up the computation in the case of the Navier-Stokes
equations is the very near of the one explained for user-defined equations
[here](@ref cs_ug_cdo_hho_base). One only focuses on the specifities arising
from the Navier-Stokes case.


Settings done in cs_user_model()
===================

The activation of the NavSto module is done thanks to the function \ref cs_navsto_system_activate
For instance, here are two examples to activate and set the main parameters for the
NavSto module

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_activate

or

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_activate2

The first parameter is the structure managing the domain boundaries. An example
of settings for the domain boundaries is described [here](@ref ug_cdo_sec_domain_boundaries).

The second parameter specifies the type of model to
consider among the following choice:
    - \ref CS_NAVSTO_MODEL_STOKES,
    - \ref CS_NAVSTO_MODEL_OSEEN,
    - \ref CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,

The *base* model can be updated thanks to the third parameter which is a flag
built from the following elemental bit (\ref cs_navsto_param_model_bit_t):
    - \ref CS_NAVSTO_MODEL_STEADY to specify the model is steady-state (by
      default, this is not the case).
    - \ref CS_NAVSTO_MODEL_GRAVITY_EFFECTS
    - \ref CS_NAVSTO_MODEL_BOUSSINESQ

The fourth parameter specifies which type of velocity-pressure algorithm will be
used (\ref cs_navsto_param_coupling_t). The choice is done among:
    - \ref CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY (cf. \cite MiBoE22 and
      \cite Mila20),
    - \ref CS_NAVSTO_COUPLING_MONOLITHIC
    - \ref CS_NAVSTO_COUPLING_PROJECTION (work in progress),

The last parameter specifies predefined post-processing operations. As the third
parameter, this a flag built from the following elemental bit
(\ref cs_navsto_param_post_bit_t):
    - \ref CS_NAVSTO_POST_VELOCITY_DIVERGENCE
    - \ref CS_NAVSTO_POST_KINETIC_ENERGY
    - \ref CS_NAVSTO_POST_VORTICITY
    - \ref CS_NAVSTO_POST_VELOCITY_GRADIENT
    - \ref CS_NAVSTO_POST_STREAM_FUNCTION (This adds an equation named
      *streamfunction_eq* and its variable field named *stream_function*. This
      equation relies on a CDO vertex-based equation.)
    - \ref CS_NAVSTO_POST_HELICITY
    - \ref CS_NAVSTO_POST_ENSTROPHY
    - \ref CS_NAVSTO_POST_MASS_DENSITY
    - \ref CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE
    - \ref CS_NAVSTO_POST_PRESSURE_GRADIENT




Predefined equations associated to the Navier--Stokes equations are
    - the *momentum* equation is automatically added when activated with a
    **monolithic** or **artificial compressibility** velocity-pressure coupling



In all cases, a vector-valued field named _velocity_ and a scalar-valued field
named _"pressure"_ are created. Moreover, several properties are added:
    - the property _mass_density_ (a \ref cs_property_t structure);
    - the property _laminar viscosity_ (a \ref cs_property_t structure);
    - the properties _turbulent_viscosity_ and the _total_viscosity_ are added
      if the model of turbulence is different from the laminar one (cf. \ref cs_turb_model_t);
    - the property _graddiv_coef_ (a \ref cs_property_t structure) when the
**artificial compressibility** is set;

along with the advection field _mass_flux_ (a \ref cs_adv_field_t structure)


Settings done in cs_user_finalize_setup()
===================

Boundary conditions
-------------------

In the case of the NavSto module, this is done as follows (One does not access
to the equation directly since it depends on the way the velocity/pressure
coupling is done).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_bc1

where the function `_vel_def` is defined as follows

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_vel_function

Definition of source terms
---------------------

Since the equation on which this source term applies, depends on the
choice of the velocity-pressure algorithm, the way to proceed varies
slightly of the way used on a user-defined equation.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_st1

where the function `_src_def` is defined as follows

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_st_function


Settings done in cs_user_parameters() {#cs_ug_cdo_navsto_set_param}
===================

The rationale is similar to the one detailed in [this section](@ref cs_ug_cdo_hho_base_user_param).

Here are listed some specificities related to the NavSto module.

Set the strategy to solve the Navier-Stokes system when a monolithic coupling is used
-------------------

When a *monolithic* velocity-pressure is set, the linear system to solve is a
_saddle-point_ problem. This class of linear systems needs specific choices of
preconditioner/solver. The default settings is not always the optimal choice in
terms of efficiency. The settings of saddle-point problems is detailed
[here](@ref cs_ug_cdo_sles_saddle)


To go further
=============

The detailed description of CDO schemes, the mathmatical analysis and numerical
results on benchmarks are available in the following PhD thesis:

* [**PhD thesis**: "Compatible Discrete Operator schemes for the unsteady incompessible Navier-Stokes equations"][Mila20] \cite Mila20


Additional publications :
* [**Article**: "Artificial compressibility methods for the incompressible Navier-Stokes equations using lowest-order face-based schemes on polytopal meshes"][MiBoE22] \cite MiBoE22
