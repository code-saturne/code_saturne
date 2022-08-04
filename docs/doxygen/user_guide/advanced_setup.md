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

<!-- ----------------------------------------------------------------------- -->

\page advanced_specific_physics Specific physical models

Use of a specific physical model
================================

Specific physical models such as dispersed phase, atmospheric flows,
gas combustion, pulverized fuel combustion, electric arcs models, and
compressible flows can be activated using the GUI (select the Calculation
features), or by using the function \ref cs_user_model of the
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
