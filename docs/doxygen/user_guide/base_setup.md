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

\page cs_ug_base_setup Base modeling setup

The documentation for this section is in the process of migration from the pdf
documentation. It is also recommended to check the
[pdf user's guide](../../../user/user.pdf) for sections which may not have been migrated yet.

[TOC]

- \subpage base_setup_param
- \subpage base_setup_initialization
- \subpage base_setup_boundary_conditions
- \subpage base_setup_physical_properties
- \subpage base_setup_source_terms
- \subpage base_setup_porosity_head_losses
- \subpage base_setup_mass_injection
- \subpage gui_user_law_editor
- \subpage base_var_modify_end_time_step

<!-- ----------------------------------------------------------------------- -->

\page base_setup_param Setup of the main parameters

Base setup should be done in the Graphical User Interface (GUI), possibly
completed by using the user functions in \ref cs_user_parameters.c.
In the GUI, the initialization is performed by filling the parameters displayed
in figures [calculation features](@ref gui_calculation_features) to
[fluid_properties](@ref gui_fluid_props).

If some options (for example 'Mobile mesh') are activated in the GUI, additional
pages may appear. In general, when activating or deactivating an option, pages
below the one currently being visited (in the setup tree) may be shown or
hidden, but pages appearing above should not be modified, so setting up
a computation by visiting GUI pages from top to bottom is recommended,
though any page may be re-visited at any time.

The sections filled for the initialization of the main parameters
are the following:

- Mesh zone definitions: volume and surface regions required by
  further definitions must be defined here, using the appropriate selection
  criteria. In many cases, zones should simply correspond to mesh groups,
  but explicitely defining them allows merging or splitting selections using
  boolean operations with group ranges and geometric criteria when needed (see
  see [mesh zone definitions](@ref gui_volume_zone_def) example figure.

- Thermophysical model options: specific physical models, ALE mobile mesh, turbulence model, thermal model and species transport (definition of the scalars and their variances);
see figures [calculation features](@ref gui_calculation_features) to  [species](@ref gui_species). If a thermal model is activated, two other sections on conjugate heat transfer and radiative transfers can be filled in (see [thermal scalar](@ref gui_thermal_scalar)).

- Body forces: gravity and coriolis forces, see [body forces](@ref gui_body_forces).

- Physical properties: reference pressure, fluid properties (density, viscosity, thermal conductivity, specific heat and scalar diffusivity), see [fluid properties](@ref gui_fluid_props).

- Volume conditions: definition of volume regions (for initialization, head losses and source terms, [user source terms](@ref sec_prg_usersourceterms) and [head losses](@ref sec_prg_headlosses), initialization of the variables (including scalars), see [figure](@ref gui_initialization).

- Boundary conditions: definition of boundary conditions.

- Time settings: number and type of time steps, and restart settings, see figure [time step](@ref gui_time_step).

- Numerical parameters: number and type of time steps, and advanced parameters for the numerical solution of the equations, see figures [global parameters](@ref gui_global_parameters) to [numerical parameters](@ref gui_numerical_parmeters).

\anchor gui_volume_zone_def
\image html gui_volume_zone_def.png "Definition of mesh zones"

\anchor gui_calculation_features
\image html gui_calculation_features.png "Calculation feature selection"

\anchor gui_turbulence_models
\image html gui_turbulence_models.png "Turbulence model selection"

\anchor gui_thermal_scalar
\image html gui_thermal_scalar.png "Thermal model selection"

\anchor gui_user_scal_def_init
\image html gui_user_scal_def_init.png "Transported species and scalars definition"

\anchor gui_body_forces
\image html gui_body_forces.png "Body forces definition"

\anchor gui_fluid_props
\image html gui_fluid_props.png "Fluid properties"

\anchor gui_initialization
\image html gui_initialization.png "Variables initialization"

\anchor gui_global_res_parameters
\image html gui_numerical_parameters.png "Global resolution parameters"

\anchor gui_numerical_parameters
\image html gui_numerical_parameters.png "Numerical parameters for the main variables"

\anchor gui_time_step
\image html gui_time_step.png "Time step settings"

For more details about the different parameters, please refer to the
[field keyword list](@ref field_keywords), [variable reference](@ref cs_var_dico),
and [examples](@ref cs_user_examples).

In addition to the GUI, user-defined functions may be used:
- \ref cs_user_model to select a given physical model.
- \ref cs_user_parameters to define most general or variable-based parameters
- \ref cs_user_time_moments to define time moments associated with the various fields.
- \ref linear_solvers to define linear solver options specific to any given system.
- \ref cs_user_finalize_setup to modify or variable-based settings, or define them for secondary variables based on options set on main parameters (such as additional variables associated with variable scalar diffusivity, turbulent Schmidt number, ...).

<!-- ----------------------------------------------------------------------- -->

\page base_setup_initialization Non-default variables initialization

The non-default variables initialization can be done using the GUI
or in the `cs_user_initialization` user-defined function.
At the calculation beginning, the variables are initialized
automatically by the code. Velocity is set to 0, scalars and properties either to zero
or to a reference value when available (such as the reference temperature or density)
and the turbulent variables are estimated from a reference velocity and characteristic length.

For <em>k</em>, in the <em>k-ε</em>, <em>R<sub>ij</sub>-ε</em>, v2f,
or <em>k-ω</em> models:
\f[
k = 1.5 \left(0.02 \textrm{ \texttt{uref}}\right)^2
\f]
and in <em>R<sub>ij</sub>-ε</em>:
\f[
R_{ij}=\frac{2}{3}k\delta_{ij}
\f]

For <em>k-ε</em> in the <em>k-ε</em>, <em>R<sub>ij</sub>-ε</em>, or <em>v2f</em>
models:
\f[
\varepsilon = k^{1.5} \frac{C_\mu}{\textrm{\texttt{almax}}}
\f]

For <em>ω</em>, in the <em>k-ω</em> model:
\f[
\omega = k^{0.5} \frac{1}{\textrm{\texttt{almax}}}
\f]

For <em>ϕ</em> and \f$ \overline{f} \f$ in the v2f models:
\f[
\left\{\begin{array}{ll}
\varphi = & \frac{2}{3} \\
\overline{f} = & 0
\end{array}\right.
\f]

For <em>α</em> in the EBRSM and BL-v2/k models:
\f[
\alpha = 1
\f]
For <em>ν<sub>t</sub></em> in the Spalart-Allmaras model:
\f[
\tilde{\nu}_t = 0.02 \sqrt{\frac{3}{2}} (\textrm{\texttt{uref}}) (\textrm{\texttt{almax}})
\f]

The \ref cs_user_initialization user-defined function allows if necessary to
initialize certain variables to values of various fields, whether those fields
represent solved variables (the most common case) or properties or even
local time step values closer to their estimated final values,
in order to obtain a faster convergence.

This function can also be used to modify values in areas not covered by
the initial mesh when restarting from a different mesh, or even defining
an interpolation from data computed with a different computational code,
as shown in the [examples](@ref user_initialization_remapper_3d).

\note **Value of the time step**

- For calculations with constant and uniform time step
  the time step is equal to the reference time step
  (\ref cs_glob_time_step->dt_ref in user functions).

- For calculations with a non-constant time step,
  the value the reference time step. In the case of a restart
  using the same time stepping options,
  the current time step should be read from the restart.

  As \ref cs_user_initialization is called after reading the restart
  file, it could be used to modify the reference or local time step
  it needed.

<!-- ----------------------------------------------------------------------- -->

\page base_setup_boundary_conditions Manage boundary conditions

As usual, except for advanced models where a high level of automation
is needed through specific preprocessing scripts or tools, using
the GUI to define boundary conditions is recommended.

Check the [pdf user's guide](../../../user/user.pdf) for details
on legacy boundary condition settings and types.

Boundary condition zones
------------------------

Prior to defining boundary conditions, the appropriate zones should be defined.
Using the GUI, this is done in the Mesh section, as an early definition
of the defined zones (before the mesh is actually read or built) may be useful
for many other settings, not limited to boundary conditions.

\anchor gui_bc_zones
\image html gui_bc_regions.png "GUI: definition of the boundary zones"

In most cases, the names of mesh groups than can be used for zone definitions
faces may be read directly from `preprocessor.log` file created by the
Preprocessor. Following an initial mesh proprocessing or verification run, these
definitions can be imported directly by the GUI under the
"Mesh/Boundary zones/Add from Preprocessor log" section.

For advanced cases, or when it is necessary to define zones which overlap,
or whose section is based on a (possibly time-varying) function,
the \ref cs_user_zones function may be used.
A few [examples](@ref cs_user_zones_examples) are also provided.

Boundary condition values
-------------------------

Based on the defined zones, the GUI can be used to define boundary
types and conditions:

\anchor gui_bc_definitions
\image html gui_bc_parameters.png "GUI: Base boundary condition definitions"

For advanced cases, user-defined functions are available, as usual.
Many [examples](@ref cs_user_boundary_conditions_examples) are provided.

## Zone-based user-defined function definitions

To define or re-define zone-based boundary condition values,
the \ref cs_user_boundary_conditions_setup (or alternatively
\ref cs_user_finalize_setup_wrapper) functions may be used.

Note that at walls and when using wall laws (which is the case with most
turbulence models), the boundary values prescribed through Dirichlet
or Robin conditions are not currently directly appplied, but used in conjunction
with the wall model. To force true Dirichlet or exchange conditions,
using the legacy boundary condition settings is still needed.

## Legacy user-defined function definitions

For definitions using the legacy system, the \ref cs_user_boundary_conditions
(in C) or \ref cs_f_user_boundary_conditions (Fortran) functions may be used.

Boundary conditions with LES
----------------------------

### Synthetic Eddy Method

The \ref cs_user_les_inflow.c user-defined function allows to generate the
unsteady boundary conditions for the LES by the Synthetic Eddy
Method.
The basic principle of this method is illustrated in the following figure.

\anchor sem_principle
\image html sem_principle.png "Synthetic Eddy Method principle"

In the figure above, \f$\mathcal{S}\f$ is the inlet boundary,
\f$\mathcal{B}\f$ the virtual box and \f$ \mathbf{U}_c \f$ the
advection velocity of the eddies.

The turbulent fluctuations at the inlet are generated
by a set of synthetic eddies advected across the inlet
boundaries. The eddies evolve in a virtual "box" surrounding the inlet
boudaries and each of them contributes to the normalized velocity
fluctuations, depending on its relative position with the inlet faces
and on a form function characterizing the shape of the eddies. By this
way, the Synthetic Eddy Method provides a coherent flow with a target
mean velocity and target Reynolds stresses at LES inlet.

\warning
LES inlets can only be defined for inlet boundary zone types.
if Dirichlet values are provided for these zones in the GUI
or user-defined functions, they are overwritten by
those provided by the Synthetic Eddy Method.

In the current version of code_saturne, the Synthetic Eddy Method is not
available through the GUI but only through the
\ref cs_user_les_inflow.c user file.

- \ref cs_user_les_inflow_define (required): define parameters of synthetic turbulence
  at LES inflow.
- \ref cs_user_les_inflow_update (advanced): update of the characteristics of a given synthetic turbulence inlet.
- \ref cs_user_les_inflow_advanced:
  definition of mean velocity, Reynolds stresses and dissipation rate
  for each boundary face of the given synthetic turbulence inlet.

Use of these functions is illustrated in the
[generation of synthetic turbulence at LES inlets](@ref les_inflow)
page.

The number of synthetic eddies in the "box" might be adjusted, depending
on the case (in particular the size of the inlet plane and the level
of turbulence). As a general rule, the greater is the better since an
insufficient number can lead to an intermittent signal while some numerical
tests have shown that this parameter does not have a great influence
beyond a threshold value. Given the inlet of size <em>h<up>2</up><em> of
a shear flow at a given Reynolds number \f$Re=u_\tau h/\nu\f$, an appropriate
number of eddies can be evaluated by \f$(Re/50)^3\f$ (<em>Re</em> and
50 approximates respectively the size, in wall unit, of the largest and
the smallest synthetic eddy.

\note Size of the synthetic eddies
The specification of the dissipation rate <em>ε</em> at the inlet is
used to compute the size \f$\sigma_i\f$ of the synthetic eddies in
the <em>i</em> Cartesian direction. One has:
\f[
\sigma_i=\max\left\{C\frac{\big(\frac{3}{2}R_{ii}\big)^{3/2}}{\varepsilon},\Delta\right\},\qquad
C=0.5.
\f]
\f$\Delta\f$ is a reference size of the grid, in order to assume that all
synthetic eddies are discretized. In the implementation of code_saturne, it is
computed at each inlet boundary face <em>F</em> as:
\f[
\Delta=2\max_{i\le3,V\in\mathcal{V}}\Big\{\big|x_i^V-x_i^C\big|\Big\}
\f]
with \f$\mathcal{V}\f$ the subset of the vertices of the boundary face <em>F</em>
and <em>C</em> the cell adjacent to <em>F</em>.

### Other LES inflow methods

For the sake of comparison, other LES inflow methods are
available, in addition to the Synthetic Eddy Method:

- The **Batten method**.\n
  With this method, the inflow velocity signal is the superposition
  of several Fourier modes. As for the Synthetic Eddy Method, the mean
  velocity, the turbulent kinetic energy and the dissipation rate have
  to be specified at the inlet: either giving their reference values
  through \ref cs_user_les_inflow_define, or their local values
  with \ref cs_user_les_inflow_advanced.

- **Random**.\n
  Turbulent fluctuations are given by a Gaussian
  noise. Only the mean velocity and Reynolds stresses need to be
  specified. The turbulent fluctuations provided by this method are
  much less realistic than those provided by the Synthetic Eddy Method
  or the Batten method. Especially for low Reynolds number flows,
  this could lead to the rapid dissipation
  of this fluctuations and the laminarization of the flow.

- **Laminar**.\n
  Adding no fluctuation, this method does not require
  any parameter. It should be reserved to regions where the flow is
  laminar.

<!-- ----------------------------------------------------------------------- -->

\page base_setup_physical_properties Manage the variable physical properties

Basic variable physical properties {#cs_ug_phys_prop}
==================================

When the fluid properties are not constant, the user may define the variation laws
in the GUI or in the \ref cs_user_physical_properties user-defined function,
which is called at each time step. In the GUI, in the item **Fluid properties**
under the **Physical properties** section , the variation laws are defined for the
fluid density, viscosity, specific heat, thermal conductivity and scalar diffusivity
as shown on the figures below.

\anchor fig_gui_fluid_props2
\image html gui_fluid_props.png "Physical properties - Fluid properties"

Variable (in space and time) properties can be defined using a
[formula editor](@ref cs_ug_meg_editor), described in a later section.

The validity of the variation laws is the user's responsibility, and
should be verified, particularly when non-linear laws are defined (for instance,
a third-degree polynomial law may produce negative density values).

\warning

- If the user wishes to impose a variable density or viscosity in
  \ref cs_user_physical_properties, a variable or viscosity must first
  be selected in the GUI, or `irovar` or `ivivar` respectively must
  have been set to 1 in \ref cs_glob_fluid_properties at an earlier
  stage (i.e. in \ref cs_user_parameters).

- In order to impose a physical property (<em>ρ</em>, <em>μ</em>, <em>λ</em>
  <em>C<sub>p</sub></em>), a reference value should be provided (except
  possibly for some physical models in which this is predefined.

- By default, the <em>C<sub>p</sub></em> coefficient and the
  diffusivity for the scalars (<em>λ</em> for the temperature) are considered
  as constant in time and uniform in space, with the \ref cs_fluid_properties_t::cp0
  and the value associated to the field's \ref diffusivity_ref keyword.
  To assign a variable value to <em>C<sub>p</sub></em>, the user **must** specify
  it in the GUI (with a user law) or set \ref cs_fluid_properties_t::icp to
  1 for \ref cs_glob_fluid_properties in \ref cs_user_parameters.c,
  and assign a value for each cell to the array `cpro_cp` which can be
  accessed through \ref CS_F_(cp)->val.

- In the same manner, to have a variable diffusivity for a given
  scalar}, the user **must** specify it in the GUI (with a user law)
  or set that field's \ref diffusivity_id keyword to a value > -1
  in \ref cs_user_parameters.c before assigning values to the matching field.

- For variable (in space) properties, it is always possible to assign a value
  in the GUI and overwrite it in \ref cs_user_physical_properties.

- The scalar \ref diffusivity_id must not be defined for user scalars
  representing the average of the square of the fluctuations of another
  scalar, because the diffusivity of a user scalar representing the average
  of the square of the fluctuations of another scalar is based on the
  diffusivity of this parent scalar.

Turbulence model properties
===========================

Modification of the turbulent viscosity
---------------------------------------

The \ref usvist user-defined subroutine can be used to modify the calculation
of the turbulent viscosity, *i.e.* <em>μ<sub>t</sub></em> in \f$kg.m^{-1}.s^{-1}\f$.
The correspondig field, `turbulent_viscosity`, can be accessed by calling
`field_get_val_s(ivisct, cpro_visct)`. The
subroutine is called at the beginning of every time step, after the
calculation of the physical parameters of the flow and of the
"conventional" value of <em>μ<sub>t</sub></em> corresponding to the chosen
turbulence model.

\warning: The calculation of the turbulent viscosity being a
particularly sensitive aspect, bad choices here may
seriously distort the results.

Modification of the variable C of the dynamic LES model
-------------------------------------------------------

The \ref cs_user_physical_properties_smagorinsky_c user-defined
function can be used to modify the calculation
of the variable *C* of the LES sub-grid scale dynamic model.

It worth recalling that the LES approach introduces the notion of
filtering between large eddies and small motions. The solved variables
are said to be filtered in an "implicit" way. Sub-grid scale models
("dynamic" models) introduce in addition an explicit filtering.

The notations used for the definition of the variable *C* used in the
dynamic models of code_saturne are specified below. These notations are
the ones assumed in the document \cite benhamadouche:tr:01.

The value of *a* filtered by the explicit filter (of width
\f$\widetilde{\overline{\Delta}}\f$) is called \f$\widetilde{a}\f$ and the value
of *a* filtered by the implicit filter (of width \f$\overline{\Delta}\f$) is
called \f$\overline{a}\f$.
We define:
\f[
\begin{array}{ll}
\overline{S}_{ij}=\frac{1}{2}(\frac{\partial \overline{u}_i}{\partial x_j}
                  +\frac{\partial \overline{u}_j}{\partial x_i})  &
||\overline{S}||=\sqrt{2 \overline{S}_{ij}\overline{S}_{ij}}\\
\alpha_{ij}=-2\widetilde{\overline{\Delta}}^2
             ||\widetilde{\overline{S}}||
             \widetilde{\overline{S}}_{ij}&
\beta_{ij}=-2\overline{\Delta}^2
             ||\overline{S}||
               \overline{S}_{ij}\\
L_{ij}=\widetilde{\overline{u}_i\overline{u}_j}-
 \widetilde{\overline{u}}_i\widetilde{\overline{u}}_j&
M_{ij}=\alpha_{ij}-\widetilde{\beta}_{ij}\\
\end{array}
\f]

In the framework of LES, the total viscosity (molecular + sub-grid) in
$kg.m^{-1}.s^{-1}$ may be written in \CS:
\f[
\begin{array}{llll}
\mu_{\text{total}}&=&\mu+\mu_{\text{sub-grid}} &
    \text{\ \ if\ \ }\mu_{\text{sub-grid}}>0\\
                   &=&\mu                          &
    \text{\ \ otherwise }\\
\text{with\ }\mu_{\text{sub-grid}}&=&\rho C \overline{\Delta}^2 ||\overline{S}||
\end{array}
\f]

\f$\overline{\Delta}\f$ is the width of the implicit filter, defined at the
cell \f$\Omega_i\f$ by\n
\f$\overline{\Delta}=XLESFL*(ALES*|\Omega_i|)^{BLES}\f$

In the case of the Smagorinsky model, *C* is a
constant which is worth \f$C_s^2\f$. \f$C_s^2\f$ is the so-called Smagorinsky
constant and is stored in the variable \ref cs_turb_csmago.

In the case of the dynamic model, *C* is variable in
time and in space. It is determined by
\f$\displaystyle C=\frac{M_{ij}L{ij}}{M_{kl}M_{kl}}\f$.

In practice, in order to increase the stability, the code does not use the
value of *C* obtained in each cell, but an average with the values
obtained in the neighboring cells (this average uses the extended
neighborhood or other vertex neighbors and corresponds to the explicit
filter). By default, the value calculated by the code is
\f[
C=\frac{\widetilde{M_{ij}L{ij}}}{\widetilde{M_{kl}M_{kl}}}
\f]

The \ref cs_user_physical_properties_smagorinsky_c function (called at each
time step, only when this model is active) allows to modify this value. It is
for example possible to compute the local average after having computed the
\f[
C=\widetilde{\left[\frac{M_{ij}L{ij}}{M_{kl}M_{kl}}\right]}
\f]

<!-- ----------------------------------------------------------------------- -->

\page base_setup_source_terms User source terms

Check [pdf user's guide](../../../user/user.pdf) for details.

<!-- ----------------------------------------------------------------------- -->

\page base_setup_porosity_head_losses Pressure drops (porosity) and head losses

Head Losses
===========

Pressure drops can be defined in the GUI or in the user sources.
In the GUI, the **Volume zones** allows to define areas where pressure drops are
applied, see an example in the [head losses region figure](@ref fig_gui_head_loss_regions).
The item **Head losses** allows to specify the head loss coefficients
(see [head loss coefficients figure](@ref fig_gui_head_loss_coeffs).
The tensor representing the pressure drops is supposed to be symmetric and positive.

\anchor fig_gui_head_loss_regions
\image html gui_head_loss_regions.png "Creation of head losses region"

\anchor fig_gui_head_loss_coeffs
\image html gui_head_loss_coeffs.png "Head loss coefficients"

In the user sources, two files can be of use: \ref cs_user_zones.c
(called at the computation start) to define volume zones and \ref cs_user_head_losses.c
(called at each iteration) to specify the values of the head losses coefficients.
As usual, volume zones defined with the GUI are available in \ref cs_user_head_losses.c.

Porosity
========

Porous zones can be defined through the GUI in the **Volume zones** page.

Alternatively, the porous model can be activated using the
\ref cs_porous_model_set_model function in \ref cs_user_parameters.
Advanced zone porosity values can be set in the
\ref cs_user_porosity function.

See \ref cs_porosity for examples.

Porous zones are defined at the beginning of the computation, before the time
loop.

<!-- ----------------------------------------------------------------------- -->

\page base_setup_mass_injection Mass volume injection

Injection of mass directly in the volume (based on mass source terms)
can be defined for selected volume zones.
The mass conservation equation is then modified as follows:
\f[
\frac{\partial \rho}{\partial t} + div(\rho\vect{u})=\Gamma
\f]

<em>Γ</em> is the mass source term expressed in
<em>kg.m</em><sup>-3</sup><em>.s</em><sup>-1</sup>.

The presence of a mass source term modifies the evolution equation of
the other variables, too. Let \f$ \varia \f$ be any solved variable apart
from the pressure (velocity component, turbulent energy, dissipation,
scalar, ...). Its evolution equation becomes:
\f[
\frac{\Delta(\rho\varia)}{\Delta t} = ... + \Gamma(\varia^{in} - \varia)
\f]

\f$ \varia^{in} \f$ is the value of \f$ \varia \f$ associated to the  mass entering
or leaving the domain. After discretization, the equation may be written:
\f[
\rho \dfrac{\varia^{(n+1)} - \varia^{(n)}}{\Delta t} = ... + \Gamma(\varia^{in} - \varia^{(n+1)})
\f]

Mass source terms can be defined using the `cs_equation_add_volume_mass_injection`
series of functions in \ref cs_user_finalize_setup.
The value assigned to the pressure variable is the mass injection rate.

For each other variable \f$ \varia \f$, there are two possibilities:

-  We can consider that the mass is added (or removed) with the
   ambient value of \f$ \varia \f$. In this case
   \f$ \varia \f$: \f$ \varia^{in} = \varia^{(n+1)} \f$ and the equation
   of \f$ \varia \f$ is not modified (so no specific definition needs to
   be added).

-  Or we can consider that the mass is added with an
   imposed value \f$ \varia^{in} \f$ (this solution is physically correct
   only when the mass is effectively added, \f$ \Gamma > 0 \f$).

For the variance, do not take into account the scalar \f$ \varia^{in} \f$
in the environment where \f$\varia \ne \varia^{in}\f$ generates a variance source.

\subpage cs_user_volume_mass_injection

Further details and examples in the linked example page above.

<!-- ----------------------------------------------------------------------- -->

\page gui_user_law_editor GUI user law editor

A formula interpreter is embedded in code_saturne, which can be used
through the GUI.In order to call the formula editor of the GUI, click on the button:
[GUI formula button](gui_formula_button.png).

This will call a popup window similar to the following one

\anchor fig_gui_density_law
\image html gui_density_law.png "Definition of a user law for the density"

The formula editor is a window with three tabs:

- User expression

  This tab is the formula editor. At the opening of the
  window only the required symbols are displayed.
  The syntax colorization shows to the user symbols which are
  required symbols, functions, or user variables.
  Each expression uses a C-like syntax and must be closed by a semicolon (`;`).
  Compared to a true C syntax, the type of local variables does not need
  to be defined, as they are assumed to de real-valued.

  Required symbols must be present in the final user law. A
  syntax checker is used when the user clicks on the OK button.

- Predefined symbols

  There are three types of symbols

  __Useful functions__

    `cos`: cosine \n
    `sin`: sine \n
    `tan`: tangent \n
    `exp`: exponential \n
    `sqrt`: square root \n
    `log`: Napierian logarithm \n
    `acos`: arc cosine \n
    `asin`: arc sine \n
    `atan(x)`: arc tangent (arc tangent of x in radians; the return value is in the range [-π/2, π/2])\n
    `atan2(y,x)`: arc tangent (arc tangent of y/x in radians; the return value is in the range [-π, π]) \n
    `cosh`: hyperbolic cosine \n
    `sinh`: hyperbolic sine \n
    `tanh`: hyperbolic tangent \n
    `abs`: absolute value \n
    `mod`: modulo \n
    `int`: floor \n
    `min`: minimum \n
    `max`: maximum \n\n

  __Useful constants__

    `pi` = 3.14159265358979323846 \n
    `e` = 2.718281828459045235\n\n

  __Operators and statements__

    `+` `-` `*` `/` `^`
    `!` `<` `>` `<=` `>=` `==` `!=` `&&` `||`
    `while` `if` `else`

- Examples

  This tab displays examples, which can be copied, pasted, and adapted.

<!-- ----------------------------------------------------------------------- -->

- \page base_var_modify_end_time_step Modification of the variables at the end of a time step

The \ref cs_user_extra_operations function is called at the end every time step.
It can be used to output of modify any variable at the end of every time step.

Several [examples](@ref cs_user_extra_operations_examples) are provided.

\warning

As all the variables (solved variables, physical
properties, geometric parameters) can be modified in this function, a
wrong use may totally distort the calculation. **If you cannot explain the
theory behind a variable value modification in this function, do not do it**.
