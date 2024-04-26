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

\page cs_ug_cdo_hho Using CDO/HHO schemes

[TOC]

<!--
    References used in this page
-->

[Bonel14]: https://hal.archives-ouvertes.fr/tel-01116527
[Cant16]: https://hal.science/tel-01419312v1
[Mila20]: https://hal.science/tel-03080530v2
[MiBoE22]: https://hal.science/hal-03215118v1


Introduction {#sec_ug_cdo_hho_intro}
============

Compatible Discrete Operator (CDO) and Hybrid High-Order (HHO) schemes are
available in code_saturne as alternative space discretizations to the legacy
Finite Volumes (FV) schemes. These are new generation of space discretizations
dedicated to polyhedral meshes.  The capabilities of these schemes in terms of
physical modelling are limited compared to the legacy FV schemes.

With **CDO schemes**, the following models/equations are available:
    - **User-defined scalar-valued or vector-valued equations** with
      - unsteady term
      - advection term (the advection-field can be user-defined or one retrieves
        the mass flux from the Navier--Stokes equations)
      - diffusion term with a homogeneous or heterogeneous +
        isotropic/orthotropic or anisotropic diffusivity tensor
      - reaction term (also called *implicit source term* in the FV parts of
        this documentation)
      - source terms
    - **steady/unsteady Stokes or Navier--Stokes equations** in a laminar flow
      regime
      - Several velocity/pressure coupling algorithms are available: monolithic,
        artifical compressibility or prediction/correction algorithm
    - **Thermal systems** with the possibility to model phase-change
      (liquid/solid)
    - **Groundwater flows**: Please refer to @ref cs_gwf_cdo for
      more details
    - **Solidification process**: Please refer to @ref cs_solidification_cdo
      for more details
    - **Wall distance** computations.
    - **Maxwell module**. This is an on-going work. Up to now, electro-static
      and magneto-static problems are available.

To get started
------------

A simple example relying on the resolution of the Laplace equation (steady
isotropic diffusion equation with Dirichlet boundary conditions) in a cube with
CDO schemes is available [here](@ref cs_cdo_laplacian). This is a a good
starting point for beginers.


Case settings for CDO/HHO schemes
=============

To set-up a CDO computation, one has to update the cs_user_parameters.c file and
edit successively the following functions
    - \ref cs_user_model,
    - \ref cs_user_parameters
    - \ref cs_user_finalize_setup

In addition, volume or boundary zones useful to set properties, equations,
etc. are defined either thanks to the GUI or by editing the function \ref
cs_user_zones for more advanced definitions from the file
`cs_user_zones.c`. This is similar to what is done with the legacy FV schemes.

More precisely, in each function above-mentioned the following settings can be
done:
    - \ref cs_user_model
      + Activation of CDO/HHO schemes,
      + Definition of the domain boundaries
      - Add new user-defined equations
      - Add new user-defined advection fields
      - Add new user-defined properties
      - Activate predefined modules (GWF, solidification, Navier--Stokes,
        thermal system or wall-distance computation)
    - \ref cs_user_parameters
      +  Definition of the time stepping strategy (if different from the one
         pre-defined in th GUI)
      +  Definition of the logging/restart frequency
      + Update of the numerical settings of each equation thanks to \ref
        cs_equation_param_set (including the way to solve the linear or the
        non-linear equations).
    - \ref cs_user_finalize_setup
      + Complete the definition of advection fields
      + Complete the definition of properties.
      + Complete the definition of an equation
        - Definition of boundary conditions
        - Definition of initial conditions
        - Defintion of source terms,
        - Activate predefined terms of a user-defined equation
          - The unsteady term by adding a link between an equation to a property
          - The advection term by adding a link to an advection field
          - The diffusion term by adding a link to a property
          - The reaction term by adding a link to a property


Settings done in cs_user_model()
===================

Activation of the CDO/HHO part
---------------------

The very first step is to activate the CDO module in the function
\ref cs_user_model There are two ways to switch on CDO/HHO schemes:
    - \ref CS_DOMAIN_CDO_MODE_ONLY for a usage of CDO or HHO in stand-lone
      (i.e. without the legacy Finite Volume approach)
    - \ref CS_DOMAIN_CDO_MODE_WITH_FV for a usage which can share some
    equations/models solved with CDO/HHO schemes and some other equations/models
    solved with the legacy Finite Volume approach.

CDO/HHO schemes can be activated within this function as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_activation

Domain boundaries {#ug_cdo_sec_domain_boundaries}
-------------------

Domain boundaries are useful for the (Navier-)Stokes equations or for the
 computation of the wall distance. Several types of domain boundaries can be
 defined. They are gathered in \ref cs_domain_boundary_type_t

The definition of the domain boundaries for CDO/HHO schemes is performed in two
steps: (1) set the default boundary and then add other boundaries which do not
fit the default one. The two possible default domain boundaries are
\ref CS_BOUNDARY_WALL or \ref CS_BOUNDARY_SYMMETRY Here is a first example.

\snippet cs_user_parameters-cdo-condif.c param_cdo_domain_boundary

Here is a second example.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_boundary


Activate predefined equations
-----------------

Several predefined models are available (and thus equation or set of
equations). When activating a predefined module, this implies:
    - adding equations and variable fields;
    - adding properties (and sometimes an advection field);
    - adapting the default numerical settings of the added equations to get more
      relevant settings.

### Navier-Stokes (NavSto module)

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



In all cases, a vector-valued field named "*velocity*" and a scalar-valued field
named "*pressure*" are created. Moreover, several properties are added:
    - the property "*mass_density*" (a \ref cs_property_t structure);
    - the property "*laminar viscosity*" (a \ref cs_property_t structure);
    - the properties "*turbulent_viscosity*" and the "*total_viscosity*" are added
      if the model of turbulence is different from the laminar one (cf. \ref cs_turb_model_t);
    - the property "*graddiv_coef*" (a \ref cs_property_t structure) when the
**artificial compressibility** is set;

along with the advection field "*mass_flux*" (a \ref cs_adv_field_t structure)


### Thermal system

The activation of the thermal module yields the following actions:
- Add an equation called "thermal_equation" (the name of the associated variable
  depends on the choice of the variable, by default, the variable is
  "temperature")
- Add the properties "mass_density", "thermal_capacity" and "thermal_conductivity"

Here is the simplest example to activate the thermal module.

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_activate_thermal_solver

The first parameter is a flag to describe the thermal model to consider. This
flag can be built from the following tags (\ref cs_thermal_model_type_bit_t)

- \ref CS_THERMAL_MODEL_STEADY
- \ref CS_THERMAL_MODEL_NAVSTO_ADVECTION

To specify the choice of the variable used in the thermal model (by default, the
temperature in Kelvin). This can be modified by adding the tag

- \ref CS_THERMAL_MODEL_IN_CELSIUS

and to change the main variable use the one of the following tag

- \ref CS_THERMAL_MODEL_USE_ENTHALPY
- \ref CS_THERMAL_MODEL_USE_TOTAL_ENERGY

The two previous options are not fully implemented.

The second parameter is a flag related to the numerical aspects. There is no tag
available up to now.

The third parameter is a flag related to the activation of automatic
post-processings.

- \ref CS_THERMAL_POST_ENTHALPY


### Groundwater flows (GWF)

Please refer to the dedicated documentation available [here](@ref cs_gwf_cdo).

### Solidification/melting module

Please refer to the dedicated documentation available [here](@ref cs_solidification_cdo).

### Wall-distance computation

It is possible to activate the computation of the distance to the wall using CDO
schemes (CDO vertex-based schemes are used) as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_wall_distance



Add a user-defined equation
-------------------------

User-defined equation with CDO/HHO schemes are added thanks to a call to the
function \ref cs_equation_add_user in \ref cs_user_model Here are several
examples:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_equation


Add a property
---------------

User-defined properties are added thanks to a call to the function
\ref cs_property_add in \ref cs_user_model Here are several examples:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_properties

The function \ref cs_property_add returns a pointer to a \ref cs_property_t
structure which can be used to set advanced parameters. If the pointer to a \ref
cs_property_t structure is not available, this can be easily recover thanks to a
call to the function \ref \cs_property_by_name

To enable the computation of the Fourier number related to a given property in
an unsteady simulation proceed as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_properties_opt


Add an advection field
---------------------

The definition of an advection field allows one to handle flows with a frozen
velocity field or the transport of scalar quantities without solving the
Navier-Stokes system. The add of a new user-defined advection field with CDO/HHO
schemes is specified as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_adv_field

When an advection is defined, it is possible to retrieve it and then set a
post-processing operation as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_adv_field_post



Settings done in cs_user_finalize_setup()
===================

Boundary conditions
-------------------

### User-defined equations

\snippet cs_user_parameters-cdo-condif.c param_cdo_setup_bcs

### NavSto module

In the case of the NavSto module, this is done as follows (One does not access
to the equation directly since it depends on the way the velocity/pressure
coupling is done).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_bc1

### Thermal system

The mechanism is the same as for user-defined equations. The name of the
equation is defined in the macro \ref CS_THERMAL_EQNAME

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_define_thermal_bc


Initial conditions
------------------

### Thermal system

The mechanism is detailed for the thermal system but the same mechanism can be
used for user-defined equations. In our example, the name of the equation is
defined in the macro \ref CS_THERMAL_EQNAME

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_define_thermal_ic

where the function `_initial_temperature` has a predefined prototype (cf. the
definition of the function pointer \ref cs_analytic_func_t)

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_initial_temperature_function


### GWF module

Here is another example extracted from \ref cs_user_parameters-cdo-gwf.c but
this is readily applicable to any equation.

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_set_ic

where the function `get_tracer_ic` has a predefined prototype (cf. the
definition of the function pointer \ref cs_analytic_func_t)

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_get_tracer_ic


Definition of properties
------------------

### User-defined properties

When a property has been added, the second step is to define this
property. According to the type of property (isotropic, orthotropic or
anisotropic) definitions differ.  Here are two examples:

\snippet cs_user_parameters-cdo-condif.c param_cdo_setup_property

### Thermal system

The mechanism is the same as for user-defined properties. The named are
predefined and associated to the following macros:

- \ref CS_THERMAL_LAMBDA_NAME for the thermal conductivity
- \ref CS_THERMAL_CP_NAME for the thermal capacity
- \ref CS_PROPERTY_MASS_DENSITY for the mass density. This property is also used
  when the NavSto module is activated.

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_define_thermal_properties


Definition of an advection field
-----------------

When an advection field has been added, the second step is to define this
advection field. Here are is an example of definition using an anlytic function
and the activation of optional features:

\snippet cs_user_parameters-cdo-condif.c param_cdo_setup_advfield


Definition of source terms
---------------------

### User-defined equation

A first simple example.

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_simple_source_terms

The second example shows an advanced usage relying on an analytic function and
an input structure. Moreover, a more accurate quadrature rule is specified.

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_source_terms

The user-defined function to compute the source term is specified as follows

\snippet cs_user_parameters-cdo-condif.c param_cdo_condif_analytic_st

and the function for the memory management of a \ref cs_xdef_t structure is

\snippet cs_user_parameters-cdo-condif.c param_cdo_condif_free_input

### NavSto module

In the case of the NavSto module, this is done as follows. Since the equation on
which this source term applies, depends on the choice of the velocity-pressure
algorithm, the way to proceed varies slightly of the way used on a user-defined
equation.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_st1

where the function *_src_def* is defined as follows

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_st_function



Add diffusion, advection, etc. to a user-defined equation
---------------------

Add terms like diffusion term,/ advection term, unsteady term, reaction terms or source terms.

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_terms






Settings done in cs_user_parameters()
===================

Change the logging frequency
----------------

The management of the level and frequency of details written by the code can be
specified for CDO/HHO schemes as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_domain_output


Change the time stepping strategy
----------------

The management of the time step with CDO/HHO schemes can be specified as
follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_time_step


Change the numerical parameters associated to an equation
-------------------------

To modify the numerical settings, the mechanism relies on a `(key, value)`
rationale.  The function to use is \ref cs_equation_param_set For instance, with
an equation called `"MyEquation"`, and a key called `CS_EQKEY_PARAM_TO_SET` with
the value `"value_to_set"`.

```
cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEquation");

cs_equation_param_set(eqp, CS_EQKEY_PARAM_TO_SET, "value_to_set");
```

If one wants to modify the settings for all equations, this is possible using
the function \ref cs_equation_set_default_param

```
cs_equation_set_default_param(CS_EQKEY_PARAM_TO_SET, "value_to_set");
```

All available keys are listed and described with \ref cs_equation_key_t

### Set the space discretization scheme

The key is \ref CS_EQKEY_SPACE_SCHEME with the possible value
- `"cdo_vb"` or `"cdovb"` to use a **CDO vertex-based discretization** (degrees
  of freedom are located at the mesh vertices). One value per vertex in the case
  of scalar-valued equation and three values per vertex in the case of a
  vector-valued equation.
- `"cdo_vcb"` or `"cdovcb"` to use a **CDO vertex+cell-based discretization**
  (degrees of freedom are located at the mesh vertices and at the mesh
  cells). One value per vertex and per cell in the case of scalar-valued
  equation and three values per vertex and pêr cell in the case of a
  vector-valued equation. Thanks to a static condensation operation, the
  algebraic system is reduced to only vertex unknows.
- `"cdo_fb"` or `"cdofb"` to use a **CDO face-based discretization** (degrees of
  freedom are located at interior and boundary faces and at mesh cell). One
  value per face and mesh cell in the case of scalar-valued equation and three
  values per face and mesh cell in the case of a vector-valued equation. Thanks
  to a static condensation operation, the algebraic system is reduced to only
  face unknows.
- `"cdo_cb"` or `"cdocb"` to use a **CDO cell-based discretization** (degrees of
  freedom are located at mesh cells for the potential and at faces for the
  flux). Only scalar-valued equation are possible. One value per cell for the
  potential. One value per face for the flux unknown (the normal component of
  the flux).
- `"cdo_eb"` or `"cdoeb"` to use **CDO edge-based discretization** (degrees of
  freedom are located at mesh edges, one scalar per edge corresponding to the
  circulation). Only vector-valued equation are handled with this
  discretization.
- `"hho_p0"` to use **HHO(k=0)** discretization relying on \f$P_0\f$ polynomial
  approximation. (degrees of freedom are located at interior and boundary faces
  and at mesh cells). One value per face and per mesh cell in the case of
  scalar-valued equation and three values per face and mesh cell in the case of
  a vector-valued equation.  Thanks to a static condensation operation, the
  algebraic system is reduced to only face unknows.
- `"hho_p1"` to use **HHO(k=1)** discretization relying on \f$P_1\f$ polynomial
  approximation. (degrees of freedom are located at interior and boundary faces
  and at mesh cells). Three values per face and four values per cell in the case
  of scalar-valued equation and nine values per face and 12 values per cell in
  the case of a vector-valued equation.  Thanks to a static condensation
  operation, the algebraic system is reduced to only face unknows.
- `"hho_p2"` to use **HHO(k=2)** discretization relying on \f$P_2\f$ polynomial
  approximation. (degrees of freedom are located at interior and boundary faces
  and at mesh cells). Six values per face and ten values per cell in the case of
  scalar-valued equation and 18 values per face and 30 values per cell in the
  case of a vector-valued equation. Thanks to a static condensation operation,
  the algebraic system is reduced to only face unknows.

\snippet cs_user_parameters-cdo-condif.c param_cdo_numerics

More details can be found in \cite Bonel14 for CDO-Vb, CDO-Fb, CDO-Cb and
CDO-Eb. CDO-Fb are also detailed in \cite Mila20. CDO-VCb are detailed in
\cite Cant16


### Set the advection scheme

Modifiy the numerical parameters related to a given user-defined equation. This
way to proceed applies to all equations: predefined or user-defined.

\snippet cs_user_parameters-cdo-condif.c param_cdo_conv_numerics



### Modifiy the numerical scheme for the diffusion term

\snippet cs_user_parameters-cdo-condif.c param_cdo_diff_numerics


### Set the linear solver and its associated options for a user-defined equation

\snippet cs_user_parameters-cdo-condif.c param_cdo_sles_settings

If the external library PETSc is available, one uses its algebraic multigrid,
otherwise one uses the in-house solver. Here an conjugate gradient (CG) with a
diagonal preconditionning (`jacobi`).

### Set the strategy to solve the Navier-Stokes system

Here is another example settings a strategy to solve a saddle-point problem
arising from a Navier-Stokes equation with a monolithic velocity-pressure
coupling. One assumes that the library is an installed dependency (see the
installation guide for more details).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_sles_alu

This is an advanced usage of the MUMPS settings.


To go further
=============

The detailed description of CDO schemes, the mathmatical analysis and numerical
results on benchmarks are available in the following PhD thesis:

* [**PhD thesis**: "Compatible Discrete Operator schemes on polyhedral meshes for elliptic and Stokes equations"][Bonel14] \cite Bonel14
* [**PhD thesis**: "Approximation of scalar and vector transport problems on polyhedral meshes"][Cant16] \cite Cant16
* [**PhD thesis**: "Compatible Discrete Operator schemes for the unsteady incompessible Navier-Stokes equations"][Mila20] \cite Mila20


Additional publications :
* [**Article**: "Artificial compressibility methods for the incompressible Navier-Stokes equations using lowest-order face-based schemes on polytopal meshes"][MiBoE22] \cite MiBoE22
