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
