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

\page cs_ug_advanced_setup Advanced modeling setup

Advanced modeling setup
=======================

The documentation for this section is in the process of migration from the pdf
documentation. It is also recommended to check the
[pdf user's guide](../../../user/user.pdf) for sections which may not have been migrated yet.

[TOC]

- \subpage advanced_specific_physics
- \subpage advanced_coal_and_gas_combution
- \subpage advanced_fuel_oil_combustion
- \subpage advanced_radiative_thermal
- \subpage advanced_conjugate_heat_transfer
- \subpage advanced_particle_tracking
- \subpage advanced_compressible
- \subpage advanced_electric_arcs
- \subpage advanced_coupling

<!-- ----------------------------------------------------------------------- -->

\page advanced_specific_physics Specific physical models

Use of a specific physical model
================================

Specific physical models such as dispersed phase, atmospheric flows,
gas combustion, pulverized fuel combustion, electric arcs models, and
compressible flows can be activated using the GUI (select the Calculation
features), or by using the \ref cs_user_model function of the
cs_user_parameters.c file (called only during the calculation initialization).

Without the GUI, the user can activate the different modules by setting the
indicators \ref cs_glob_physical_model_flag in the cs_user_model function, (see
[Base model related options examples](@ref cs_user_parameters_h_cs_user_model)
for some examples)

*WARNING: Only one specific physical module can be activated at the same
time.*

In the framework of the gas combustion modeling, users may impose
their own enthalpy-temperature tabulation (conversion law). The `indjon`
indicator must be set to 0 in this case (the default value being 1) and the name
of the tabulation file defined through the **ficfpp** variable (see
[Specific physical model activation (usppmo) examples](@ref cs_f_user_parameters_h_usppmo)).
For more details, the user may refer to the following note
(thermochemical files).

Thermo-chemical data files
--------------------------

The user must not forget to place in the `DATA` directory  the
thermochemical file `dp_C3P` (**indjon = 1**), `dp_C3PSJ` or `dp_ELE`
(depending on the selected module) Such files can be copied to the
`DATA` directory of a study case using the GUI (see the Tools menu).
Their content is described below.

For example, using gas combustion if (indjon=1) we can directly copy the
`dp_C3P` file from the GUI.

\anchor gui_dp_C3P_gas_combustion
\image html gui_dp_C3P_gas_combustion.png "Add tabulation file"

- Example of file for the gas combustion: \n
if the enthalpy-temperature conversion data base JANAF is used: `dp_C3P`
(see following array).

| Lines | Examples of values     | Variables        | Observations                                           |
|-------|------------------------|------------------|--------------------------------------------------------|
| 1     | 5                      | ngaze            | Number of current species                              |
| 2     | 10                     | npo              | Number of points for theenthalpy-temperature table     |
| 3     | 300.                   | tmin             | Lower temperature limit for the table                  |
| 4     | 3000.                  | tmax             | Upper temperature limit for the tabulation             |
| 5     |                        |                  | Empty line                                             |
| 6     | CH4 O2 CO2 H2O N2      | nomcoe(ngaze)    | List of the current species                            |
| 7     | .35 .35 .35 .35 .35    | kabse(ngaze)     | Absorption coefficient of  the current species         |
| 8     | 4                      | nato             | Number of elemental species                                                            |
| 9 \n 10 \n 11 \n 12 | .012  1  0  1  0  0 \n .001  4  0  0  2  0 \n .016  0  2  2  1  0\n .014  0  0  0  0  2 | wmolat(nato), \n atgaze(ngaze,nato) | Molar mass of the elemental species (first column) \n Composition of the current species \n as a function of the elemental species \n (ngaze following columns) |
|13     |         3              | ngazg   | Number of global species. \n Here, ngazg = 3 (Fuel, Oxidiser and Products)|
|14 \n 15 \n 16    |  1. 0. 0. 0. 0. \n 0. 1. 0. 0. 3.76 \n 0. 0. 1. 2. 7.52 | compog(ngaze,ngazg) | Composition of the global species as a \n function of the current species of line 6  In the order: Fuel (line 15), Oxidiser (line 16) and Product (line 17)  |
|17     |         1              |  nrgaz   | Number of global reactions. \n Here nrgaz = 1 (always equal to 1in this version) |
|18     | 1 2 -1 -9.52 10.52| igfuel(nrgaz), \n igoxy(nrgaz),\n stoeg(ngazg,nrgaz) | Numbers of the global species concerned by \n the stoichiometric ratio \n (first 2 integers) \n Stoichiometry in global species reaction. Negative for the reactants (here “Fuel” and “Oxidiser”) and positive for the products (here “Products”) |

if the user provides an enthalpy-temperature tabulation (there must
be three chemical species and only one reaction): `dp_C3PSJ` (see following array).
This file replaces `dp_C3P`.


| Lines | Examples of values          | Variables                       | Observations              |
|-------|-----------------------------|---------------------------------|---------------------------|
|    1  |             6                                                                       |   npo                           | Number of tabulation points |
|2 \n 3 \n 4 \n 5 \n 6 \n 7 |50. -0.32E+07 -0.22E+06 -0.13E+08 \n 250. -0.68E+06 -0.44E+05 -0.13E+08 \n 450. 0.21E+07 0.14E+06 \n 650. 0.50E+07 0.33E+06 -0.12E+08 \n 850. 0.80E+07 0.54E+06 -0.12E+08 \n 1050. 0.11E+08 0.76E+06 -0.11E+08 |th(npo), \n ehgazg(1,npo),\n ehgazg(2,npo),\n ehgazg(3,npo)| Temperature(first column), \n mass enthalpies of fuel, oxidiser \n and products (columns 2,3 and 4) \n from line 2 to line npo+1 |
|8   | .00219 .1387 .159 | wmolg(1), \n wmolg(2), \n wmolg(3) | Molar masses of fuel, \n oxidiser \n and products|
|9   | .11111   | fs(1) | Mixing rate at the stoichiometry (relating to Fuel and Oxidiser) |
|10  | 0.4 0.5 0.87 | ckabsg(1),\n ckabsg(2), \n ckabsg(3) | Absorption coefficients of the fuel, oxidiser and products |
|11  | 1. 2. | xco2, xh2o | Molar coefficients of __CO2__ \n and __H2O__ in the products \n (using Modak radiation) |



- Example of file for the electric arcs:\n
`dp_ELE` file (see following array).

| Lines | Examples of values           | Variables    | Observations                     |
|-------|------------------------------|--------------|----------------------------------|
| 1     | # Free format ASCII file ... |              | Free comment                     |
| 2     | # Comment lines ...          |              | Free comment                     |
| 3     | #...                         |              | Free comment                     |
| 4     | # Argon properties ...       |              | Free comment                     |
| 5     | # ...                        |              | Free comment                     |
| 6     | # No of NGAZG and No ...     |              | Free comment                     |
| 7     | # NGAZG NPO ...              |              | Free comment                     |
| 8     | 1 238                        | ngazg \n npo | Number of species \n Number of given temperature points for the tabulated physical properties (npo 6 npot set in ppthch) \n So there will be ngazg blocks of npo lines each |
| 9     | # ...                        |              | Free comment                     |
| 14    | 0                            | ixkabe       | Radiation options for xkabe      |
| 15    | # ...                        |              | Free comment                     |
| 16    | # Propreties ...             |              | Free comment                     |
| 17    | # T      H ...               |              | Free comment                     |
| 18    | # Temperature Enthalpy ...   |              | Free comment                     |
| 19    | # ...                        |              | Free comment                     |
| 20    | # K      J/kg ...            |              | Free comment                     |
| 21    | # ...                        |              | Free comment                     |
| 22    | # 300.    14000.  ...        | h \n roel \n cpel \n sigel \n visel \n xlabel \n xkabel | In line tabulation of the physical properties as a function of the temperature in Kelvin \n for each of the ngazg species \n Enthalpy in J/kg \n Density in kg/m3 \n Specific heat in J/(kg K) \n Electric conductivity in Ohm/m \n Dynamic viscosity in kg/(m s) \n Thermal conductivity in W/(m K) \n Absorption coefficient (radiation) |

<!-- ----------------------------------------------------------------------- -->

\page advanced_coal_and_gas_combution Pulverized coal and gas combustion module

Initialization of the variables
===============================

For **Reactive flows (combustion)**, it is possible to initialize the specific variables in the **Graphical User Interface (GUI)** or in the \ref cs_user_initialization function.

In the GUI, when a **Reactive flows (combustion)** is selected in the item **“Calculation features”**, an additional item appears: __“Gas combustion”__ the user can change it by __“Pulverized coal”__. In this item the user can define coal types, their composition, the oxidant and reactions parameters, see the following figure.

\anchor gui_coal_classes
\image html gui_coal_classes.png "Thermophysical models - Pulverized coal, coal classes"

\anchor gui_coal_composition
\image html gui_coal_composition.png "Pulverized coal combustion, coal composition"

\anchor gui_coal_reaction
\image html gui_coal_reaction.png "Pulverized coal combustion, reaction parameters"

\anchor gui_coal_oxydant
\image html gui_coal_oxydant.png "Pulverized coal combustion, oxydant"

If the user activates **gas combustion**  and does not want to use the GUI, she or he can directly use
the \ref cs_user_initialization  function (for some examples see [Initialization examples](@ref user_initialization_base_s_init)). \n

Boundary conditions
===================

**For pulverized coal**, it is possible to manage the boundary conditions in
the Graphical User Interface (GUI). When the **boundary zones** is actived,
the user specific boundary conditions are activated for inlets

\anchor gui_coal_bc
\image html gui_coal_bc.png "Boundary conditions for the pulverized coal module"


**For gas combustion** it is also possible to manage the boundary conditions for the inlet in the Graphical User Interface (GUI). The user can choose between the **burned gas** or the **Unburned gas**,
impose the the mass flow and velocity

\anchor gui_gas_bc
\image html gui_gas_bc.png "Boundary conditions for gas combustion"


Initialization of the options of the variables
==============================================

In the case of **gas combustion** or **pulverized coal combustion**, time averages, chronological records and log follow-ups can be set in GUI or in the \ref cs_user_parameters function. In the GUI, under the heading **Calculation control**, additional variables appear in the list in the items **Time averages** and **Profiles**, as well as in the item **Volume solution control**

\anchor gui_coal_time_average
\image html gui_coal_time_average.png "Calculation control - Time averages"

\anchor gui_coal_solution_control
\image html gui_coal_solution_control.png "Calculation control - Volume solution control"

For **gas combustion**, if the GUI is not used for **coal combustion**, the
\ref cs_user_parameters and \ref cs_user_combustion functions, called at calculation
start, can be used to:

- set the relaxation coefficient of the density \ref srrom.

- set the dynamic viscosity \ref diftl0.

- set the value of the constant \ref cebu of the Eddy Break Up model (only in
\ref cs_user_combustion).

<!-- ----------------------------------------------------------------------- -->

\page advanced_fuel_oil_combustion Heavy fuel oil combustion module

Initialization of transported variables
=======================================

To initialize or modify (in case of a restart) values of transported
variables and of the time step, the standard \ref cs_user_initialization function is used.

Boundary conditions
===================

Boundary conditions are defined as usual on a per-face basis in \ref cs_user_boundary_conditions or with the **GUI**.

<!-- ----------------------------------------------------------------------- -->

\page advanced_radiative_thermal Radiative thermal transfers in semi-transparent gray media

Initialization of the radiation main parameters
===============================================

The main radiation parameters can be initialized in the GUI or in the \ref cs_user_radiative_transfer_parameters user function (see [Initialization examples](@ref cs_user_radiative_transfer_h_cs_user_radiative_transfer_parameters)).
In the GUI, under the heading **Thermal models**, when one of the two thermal radiative transfers models is selected, see [Figure 1](@ref gui_rad_transf_do_params)
additional items appear. The user is asked to choose the number of directions for angular discretization, to define the absorption coefficient and specify if the radiative calculation is restarted from a checkpoint or not;
see [Figure 1](@ref gui_rad_transf_do_params) and [Figure 3](@ref gui_rad_transf_p1_params).  When **Advanced options** is selected for both models [Figure 2](@ref gui_rad_transf_do_advanced) and [Figure 4](@ref gui_rad_transf_p1_advanced)  appear, the user must fill the resolution frequency and verbosity levels. In addition, the activation of the radiative transfer leads to the creation of a **Surface solution control** item under the heading **Calculation control**, see [Figure 5](@ref gui_rad_transf_post_output), where radiative transfer variables can be selected to appear in the output log.

\warning when a calculation is run using a specific physics module,
this first heading must not be completed. The radiation module is then
activated or not, according to the parameter file related to the considered
specific physics.

Radiative transfer boundary conditions
======================================

These conditions can be defined with the GUI or by programming the \ref cs_user_radiative_transfer_bcs  function; see [boundary examples](@ref cs_user_radiative_transfer_h_boundary_conditions).
In the GUI, when one of the **Radiative transfers** options is selected in [Figure 1](@ref gui_rad_transf_do_params) and **Boundary zones** is defined in **Mesh**,
it activates specific boundary conditions each time a **Wall** is defined, see [Figure 6](@ref gui_rad_transf_wall_model). The user can then choose between 3 cases. The parameters that must be specified are displayed for one of them in [Figure 7](@ref gui_rad_transf_wall_params).

Absorption coefficient of the medium, boundary conditions for the luminance, and calculation of the net radiative flux
======================================================================================================================

When the absorption coefficient is not constant, the \ref cs_user_rad_transfer_absorption functionis called at each time
step. It is composed of three parts. In the first one, the user must provide the absorption coefficient of the medium in the array CK,
for each cell of the fluid mesh. By default, the absorption coefficient of the medium is 0, which corresponds to a transparent medium.
For more detail see [Absorption](@ref abso_flux)

\warning
When a specific physical model is activated, it is forbidden to
set the absorption coefficient in this function. In this
case, the coefficient is either calculated automatically, or provided by the user via a thermo-chemical parameter file (dp_C3P or dp_C3PSJ for gas combustion,
and dp_FCP for pulverized coal combustion).

\anchor gui_rad_transf_do_params
\image html gui_rad_transf_do_params.png "Radiative transfers - parameters of the DO method"

\anchor gui_rad_transf_do_advanced
\image html gui_rad_transf_do_advanced.png "Radiative transfers - advanced parameters of the DO method"

\anchor gui_rad_transf_p1_params
\image html gui_rad_transf_p1_params.png "Radiative transfers - parameters of the P-1 model"

\anchor gui_rad_transf_p1_advanced
\image html gui_rad_transf_p1_advanced.png "Radiative transfers - advanced parameters of the P-1 model"

\anchor gui_rad_transf_post_output
\image html gui_rad_transf_post_output.png "Calculation control - Radiative transfers post-processing output"

\anchor gui_rad_transf_wall_model
\image html gui_rad_transf_wall_model.png "Boundary conditions - choice of wall thermal radiative transfers"

\anchor gui_rad_transf_wall_params
\image html gui_rad_transf_wall_params.png "Boundary conditions - example of wall thermal radiative transfer"

<!-- ----------------------------------------------------------------------- -->

\page advanced_conjugate_heat_transfer Conjugate heat transfer

Thermal module in a 1D wall
===========================

The \ref cs_user_1d_wall_thermal function takes into account the wall-affected thermal inertia.
Some boundary faces are treated as a solid wall with a given thickness, on which the code resolves
a one-dimensional equation for the heat conduction. The coupling between the 1D module and the fluid
works in a similar way to the coupling with the **SYRTHES**. By construction, the user is not able to account
for the heat transfer between different parts of the wall. A physical analysis of each problem, case by case
is required in order to evaluate the relevance of its usage by way of a report of the simple conditions
(temperature, zero-flux ) or a coupling with **SYRTHES**.

The use of this code requires that the thermal scalar is
defined as (\ref cs_thermal_model_field() \f$ \ne \f$  NULL).

\warning
The 1D thermal module is developed assuming the thermal scalar
as a temperature. If the thermal scalar is an enthalpy, the code calls the
enthalpy to temperature conversion as defined by the model defaults,
or by the user in \texttt{cs\_user\_physical\_properties} for each
transfer of data between the fluid and the wall in order to convert the
enthalpy to temperature and vice-versa. If the thermal
variable is the total (compressible) energy, the thermal module will not work.

Internal Fluid-Thermal coupling
===============================

When at least one volume zone is defined as being solid
(see [Figure 1](@ref gui_internal_coupling_vol_zone)), scalar variables (especially
thermal scalar variables) may be solved in a fully coupled manner across the fluid
and solid domains.

For this purpose, the **Internal coupling** should be activated for the desired variables
in the matching tab of the **Coupling parameters** page, as shown in figure
[figure 2](@ref gui_internal_coupling)). This section should appear when
at least one volume zone is defined as solid.

\anchor gui_internal_coupling_vol_zone
\image html gui_internal_coupling_vol_zone.png "Solid volume zone definition"

\anchor gui_internal_coupling
\image html gui_internal_coupling.png "Conjugate heat transfer: internal coupling"

Fluid-Thermal coupling with SYRTHES
===================================

Coupling **code_saturne** with **syrthes** for **conjugate heat transfer** can be defined through
the GUI or the \ref cs_user_syrthes_coupling user function see [syrthes coupling examples](@ref cs_user_coupling_h_cs_user_syrthes_coupling).

To set such a coupling in the GUI, a thermal scalar must be
selected first in the item **Thermal scalar** under the heading **Thermophysical models**.
At least one wall boundary condition must be set to **SYRTHES coupling** type, and
the name of the associated **syrthes** instance (i.e. base directory name of the associated
solid case definition) be set, as shown in, [Figure 3](@ref gui_syrthes_coupling_bc).
The **Syrthes coupling** tab will then be available in the **Coupling parameters**
section (see [Figure 4](@ref gui_syrthes_coupling)), fo further advanced or global settings.
The zones where the coupling occurs must be defined and a projection axis can be
specified in case of 2D coupling.

\anchor gui_syrthes_coupling_bc
\image html gui_syrthes_coupling_bc.png "Boundary conditions - coupling with syrthes"

\anchor gui_syrthes_coupling
\image html gui_syrthes_coupling.png "Coupling parameters - coupling with syrthes"

<!-- ----------------------------------------------------------------------- -->

\page advanced_particle_tracking Particle-tracking (Lagrangian) module

General information
===================

 - The particle-tracking (or Lagrangian) module enables the simulation of poly-dispersed particulate flows,
   by calculating the trajectories of individual particles, mainly characterized by their diameter and density
   (if no heat nor mass transfer between particle and fluid are activated).

 - The standard use of the particle-tracking module follows the **Moments/PDF approach**: the instantaneous
   properties of the underlying flow needed to calculate the particle motion are reconstructed from the
   averaged values (obtained by Reynolds-Averaged Navier-Stokes simulation) by using stochastic processes.
   The statistics of interest are then obtained through Monte-Carlo simulation.

 - As a consequence, is is important to emphasize that the most important (and physically meaningful) results
   of a particle-tracking calculation following the Moments/PDF approach are **statistics**.
   Volume and surface statistics, steady or unsteady, can be calculated. Individual particle trajectories
   (as 1D, _EnSight_-readable cases) and displacements (as _EnSight_-readable animations) can also be provided, but only for illustrative purposes.

Activating the particle-tracking module
=======================================

The activation of the particle-tracking module is performed either:
    - in the Graphical User Interface (GUI): _Calculation features_ --> _Homogeneous Eulerian - VoF model_ --> _particles and droplets tracking_
    - or in the user function \ref cs_user_lagr_model.

Basic guidelines for standard simulations
=========================================

Except for cases in which the flow conditions depend on time, it is generally recommended to perform a first Lagrangian calculation whose aim is to reach a steady-state (i.e. to reach a time starting from which the relevant statistics do not depend on time anymore). In a second step, a calculation restart is done to calculate the statistics. When the single-phase flow is steady and the particle volume fraction is low enough to neglect the particles influence on the continuous phase behaviour, it is recommended to perform a Lagrangian calculation on a frozen field.

It is then possible to calculate steady-state volumetric statistics and to give a statistical weight higher than 1 to the particles, in order to reduce the number of simulated (**numerical**) particles to treat while keeping the right concentrations. Otherwise, when the continuous phase flow is steady, but the two-coupling coupling must be taken into consideration, it is still possible to activate steady statistics.
When the continuous phase flow is unsteady, it is no longer possible to use steady statistics. To have correct statistics at every moment in the whole calculation domain, it is imperative to have an established particle seeding and it is recommended (when it is possible) not to impose statistical weights different from the unity.

Finally, when the so-called complete model is used for turbulent dispersion modelling, the user must make sure that the volumetric statistics are directly used for the calculation of the locally undisturbed fluid flow field.

When the thermal evolution of the particles is activated, the associated particulate scalars are always the inclusion temperature and the locally undisturbed fluid flow temperature expressed in degrees Celsius, whatever the thermal scalar associated with the continuous phase is (i.e. temperature or enthalpy). If the thermal scalar associated with the continuous phase is the temperature in Kelvin, the unit is converted automatically into Celsius. If the thermal scalar associated with the continuous phase is the enthalpy, a _temperature_ property or postprocessing
field must be defined. In all cases, the thermal backward coupling of the dispersed phase on the continuous phase is adapted to the thermal scalar transported by the fluid.

Prescribing the main modelling parameters
=========================================

Use of the GUI
--------------

In the GUI, the selection of the Lagrangian module activates the heading _Particle and droplets tracking_ in the tree menu. The initialization is performed in the three items included in this heading:
    -  _Global settings_. The user defines in this item the kind of Euler/Lagrange multi-phase treatment, the main parameters, and the specific physics associated with the particles, see [Figure 1](@ref gui_lagr_global_settings)
    - _Statistics_. The user can select the volume and boundary statistics to be post-processed see [Figure 2](@ref gui_lagr_statistics).
    - _Output_. An additional entry in the postprocessing section allows defining the output frequency and post-processing options for particles and selecting the variables that will appear in the log see [Figure 3](@ref gui_lagr_output).

\anchor gui_lagr_global_settings
\image html gui_lagr_global_settings.png "Lagrangian module - View of the _Global Settings_ page"

\anchor gui_lagr_statistics
\image html gui_lagr_statistics.png "Lagrangian module - statistics"

\anchor gui_lagr_output
\image html gui_lagr_output.png "Lagrangian module - output"

Use of the function cs_user_lagr_model
--------------------------------------

When the GUI is not used, \ref cs_user_lagr_model must be completed. This function
gathers in different headings all the keywords which are
necessary to configure the Lagrangian module. The different headings refer to:
    - the global configuration parameters
    - the specific physical models describing the particle behaviour
    - the backward coupling (influence of the dispersed phase on the
      continuous phase)
    - the numerical parameters
    - the volumetric statistics
    - the boundary statistics

For more details about the different parameters and some examples, the user may refer to [examples](@ref cs_user_lagr_module_intro)

Prescribing particle boundary conditions
========================================

In the framework of the multiphase Lagrangian modelling, the management of the boundary conditions concerns the particle behaviour when there is an interaction between its trajectory and a boundary face. These boundary conditions may be imposed independently of those concerning the Eulerian fluid phase (but they are of course generally consistent). The boundary condition zones are actually redefined by the Lagrangian module ([boundary zones](@ref cs_user_lagr_boundary_conditions_h_zones)), and a type of particle behaviour is associated with each one. The boundary conditions related to particles can be defined in the Graphical User Interface (GUI) or in the \ref cs_user_lagr_boundary_conditions.c} file. More advanced user-defined boundary conditions can be prescribed in the \ref cs_user_lagr_in function from \ref cs_user_lagr_particle.c}.

Use of the GUI
--------------

In the GUI, selecting the Lagrangian module in the activates the item _Particle boundary conditions_ under the heading _Boundary conditions_ in the tree menu. Different options are available depending on the type of standard boundary conditions selected (wall, inlet/outlet, etc...),
see [Figure 3](@ref gui_lagr_bc).

\anchor gui_lagr_bc
\image html gui_lagr_bc.png "Lagrangian module - boundary conditions"

Advanced particle-tracking set-up
=================================

In this section, some information is provided for a more advanced numerical set-up of a particle-tracking simulation.

User-defined stochastic differential equations
----------------------------------------------

An adaptation in the \ref cs_user_lagr_sde function is required if
supplementary user variables are added to the particle state vector for more explanation see \ref cs_user_lagr_sde_page.

If necessary, the thermal characteristic time \f$\tau_c\f$, whose calculation can be modified by the user in the function
\ref cs_user_lagr_rt.

User-defined particle relaxation time
-------------------------------------

The particle relaxation time may be modified in the \ref cs_user_lagr_rt function according to the chosen formulation of the drag coefficient. The particle relaxation time, modified or not by the user, is available in the array _taup_ see \ref cs_user_lagr_module_time_relaxation for examples

User-defined particle thermal characteristic time
-------------------------------------------------

The particle thermal characteristic time may be modified in the \ref cs_user_lagr_rt_t function according to the chosen correlation for the calculation of the
Nusselt number see \ref cs_user_lagr_module_thermal_relaxation for examples.

<!-- ----------------------------------------------------------------------- -->

\page advanced_compressible Compressible module

When the **compressible module** is activated, it is recommended to:
    - use the option _time step variable in time and uniform in space_ (idtvar=1) with a maximum
      Courant number of 0.4 (\ref coumax = 0.4): these choices must be written in \ref cs_user_parameters.c
      or specified with the **GUI**
    - keep the convective numerical schemes proposed by default _i.e._: upwind scheme

With the compressible algorithm, the specific total energy is a new solved variable
CS_F_(e_tot). The temperature variable deduced from the specific total energy variable is
CS_F_(t_kelvin) for the compressible module.\n
Initialization of the options of the variables, boundary conditions, initialisation of the variables and
management of variable physical properties can be done with the **GUI**. We describe below the functions
the user has to fill in without the **GUI**.

Initialization of the options of the variables
==============================================

When the GUI is not being used, the function \ref cs_user_parameters in \ref cs_user_parameters.c
must be completed by the user.\n This function allows to activate the compressible (see \ref cs_user_parameters_h_cs_user_model)
module and to specify the molecular viscosity (ivivar see \ref cs_user_parameters_h_param_fluid_properties),

Management of the boundary conditions
=====================================

When running the compressible module without a GUI, the \ref cs_user_boundary_conditions function can be used to define specific boundary conditions
(see the \ref advanced_loc_var_ce file for examples of boundary conditions with the compressible module).

With the compressible module, the following types of boundary condition are avaliable:\n

  - Inlet/outlet for which velocity and two thermodynamics variables are known see \ref compressible_ex_1.
  - Supersonic output see \ref compressible_ex_2.
  - Subsonic  input with density and velocity see \ref compressible_ex_3.
  - Subsonic outlet \ref compressible_ex_4.
  - Wall (adiabatic or not) \ref compressible_ex_5.


Initialization of the variables
===============================

When the **GUI** is not used, the function \ref cs_user_initialization is used
to initialize the velocity, turbulence and passive scalars (see
the \ref user_initialization_compressible for examples of initialisations with
the compressible module). Concerning pressure, density, temperature and specific total energy, only 2 variables out
of these 4 are independent. The user may then initialise the desired variable pair
(apart from temperature-energy) and the two other variables will be
calculated automatically by giving the right value to the variable
ithvar see \ref user_initialization_comp_s_init for example.

Management of variable physical properties
==========================================

Without the **GUI**, all of the laws governing the physical properties of the fluid
(molecular viscosity, molecular volumetric viscosity, molecular thermal conductivity and
molecular diffusivity of the user-defined scalars) can be specified in the function \ref cs_user_physical_properties of
the \ref cs_user_physical_properties.c file.

The user should check that the defined laws are valid for
the whole variation range of the variables. Moreover, as only the perfect gas with a constant
adiabatic coefficient equation of state is available, it is not advised to give a law for the isobaric
specific heat without modifying the equation of state in the function \ref cs_cf_thermo which is not
a user function.

For some examples we can see:
 - [Ex. 1: molecular viscosity varying with temperature](@ref example1_comp)
 - [Ex. 2: molecular volumetric viscosity varying with temperature](@ref example2_comp)
 - [Ex. 3: isobaric specific heat varying with temperature](@ref example3_comp)
 - [Ex. 4: molecular thermal conductivity varying with temperature](@ref example4_comp)
 - [Ex. 5: molecular diffusivity of user-defined scalars varying with temperature](@ref example5_comp)


<!-- ----------------------------------------------------------------------- -->

\page advanced_electric_arcs Management of the electric arcs module

Activating the electric arcs module
===================================

The electric arcs module is activated either:

 - in the Graphical User Interface _GUI_: __Calculation features__ --> __Electrical arcs__, the user can choose between _Joule Effect_ for _joule model_ and _Joule Effect_ and _Laplace Forces_ for electric arc
 - or in the user function \ref cs_user_model in cs_user_parameters.c file, by setting the \ref cs_glob_physical_model_flag[\ref CS_ELECTRIC_ARCS] or \ref cs_glob_physical_model_flag[\ref CS_JOULE_EFFECT] parameter to a non-null value.

Initialisation of the variables
===============================

The function \re cs_user_initialization allows the user to initialise some of the specific physics variables prompted via \ref cs_user_model. It is called only during the initialisation of the calculation. As usual,the user has access to many geometric variables so that the zones can be treated separately if needed (see [Electric arcs example](@ref user_initialization_electric_arcs)).

The values of potential and its constituents are initialised if required.

It should be noted that the enthalpy is relevant.

 - For the _electric arcs_ module, the _enthalpy_ value is taken from the temperature
 of reference \ref t0 (given in \ref cs_user_parameters.c)
 from the temperature-enthalpy tables supplied in the data file **dp_ELE**.
 The user must not intervene here.

- For the _Joule effect_ module, the value of enthalpy must be specified by the user.
 Examples of temperature to enthalpy conversion are given in
 \ref cs_user_physical_properties.c). If not defined, a simple default
 law is used (\f$H = C_p T\f$).

Variable physical properties
============================

All the laws of the variation of physical data of the fluid are written (when necessary)
in the function \ref cs_user_physical_properties.

\warning
 For the _electric module_, it is here that all the physical variables are defined
 (including the relative cells and the eventual user scalars): \ref cs_user_physical_properties _is not used_.

The user should ensure that the defined variation laws are valid for the whole range of
variables. Particular care should be taken with non-linear laws (for example, a
 \f$3^{rd}\f$ degree polynomial law giving negative values of density)

\warning
 In the _electric module_, all of the physical properties are considered as variables
 and are therefore stored using the  cs_field API. \ref cp0, \ref viscls0 and \ref viscl0
 are not used

For the Joule effect, the user is required to supply the physical properties in the
function. Examples are given which are to be adapted by the user. If the temperature is
to be determined to calculate the physical properties, the solved variable, enthalpy must
 be deduced. The preferred temperature-enthalpy law should be defined
 (a general example is provided in (\ref cs_user_physical_properties),
 and can be used for the initialisation of the variables in
 (\ref cs_user_initialization)).
 For the _electric arcs_ module, the physical properties are interpolated from the data file
 __dp_ELE__ supplied by the user. Modifications are generally not necessary.

Boundary conditions
===================

Boundary conditions can be handled in the GUI or in the cs_user_boundary_conditions function as usual (see [Electric example](@ref electric_arcs_examples).
 In the \ref cs_user_boundary_conditions report, the main change from the users point of view concerns the
 specification of the boundary conditions of the potential, which isn't
 implied by default. The Dirichlet and Neumann conditions must be imposed
 explicitly using \ref icodcl and \ref rcodcl (as would be done for the classical scalar).

Furthermore, if one wishes to slow down the power dissipation (Joule
effect module) or the current (electric arcs module) from the imposed values,
they can be changed by the potential scalar as shown below:

 - For the electric arcs, the imposed current intensity can be a fixed variable and initialize by the GUI see [Figure 1](@ref gui_electric_arcs)

\anchor gui_electric_arcs
\image html gui_electric_arcs.png "Imposed current intensity"

 - For the *Joule model*, the imposed power can be a fixed variable in the same way as the electric arcs.

\warning :
*In the case of alternating current, attention should be paid to the values of potential
 imposed at the limits: the variable named "real potential" represents an affective
 value if the current is in single phase, and a "real part" if not.*

 - For the Joule studies, a complex potential is sometimes needed
 (*in the GUI Electrical model -> three-phase*): this is the  case in particular where the current
 has three phases. To have access to the phase of the potential, and not just to its
 amplitude.

 - For the Joule studies in which one does not have access to the phases, the real
 potential (imaginary part =0) will suffice (*in the GUI Electrical model -> AC/DC*): this is
 obviously the case with
 continuous current, but also with single phase alternative current. In *code_saturne*
there is only 1 variable for the potential,  called "real potential". Pay attention to
 the fact that in alternate current, the "real potential" represents a effective value
 of potential , \f$\frac{1}{\sqrt{2}}\,Pp_{max}\f$ (in continuous current there is no
 such ambiguity).

Additions for transformers
--------------------------

The following additional boundary conditions must be defined for tansformers:
  - the intensity at each electrode
  - the voltage on each terminal of transformers. To achieve it, the intensity,
    the rvoltage at each termin, the Rvoltage, and the total intensity of the
    transformer are calculated.

Finally, a test is performed to check if the offset is zero or if a boundary
 face is in contact with the ground.

<!-- ----------------------------------------------------------------------- -->

\page advanced_coupling coupling with saturne

code saturne-code saturne coupling
==================================

The user function \ref cs_user_saturne_coupling in \ref cs_user_coupling.c is
used to couple *code_saturne* with itself.
It is used for *turbo-machine* applications for instance, the first *code_saturne* managing
the fluid around the rotor and the other the fluid around the stator.
In the case of a coupling between two *code_saturne* instances, first argument *saturne_name*
of the function \ref cs_sat_coupling_define is ignored.
In case of multiple couplings, a coupling will be matched with available *code_saturne*
instances based on that argument, which should match the directory name for the
given coupled domain see [examples](@ref cs_user_coupling_h_cs_user_saturne_coupling).

Fluid-Structure external coupling
=================================

The function \ref usaste belongs to the module dedicated to external
Fluid-Structure coupling with *Code_Aster*. Here one defines the boundary
faces coupled with *Code_Aster* and the fluid forces components which are
given to structural calculation. When using external coupling with *Code_Aster*,
structure numbers necessarily need to be negative;\n
the references of coupled faces being (*i.e. -1, -2*), etc.
For examples on the function we can see [examples](@ref cs_user_fluid_structure_interaction_h_usaste)
