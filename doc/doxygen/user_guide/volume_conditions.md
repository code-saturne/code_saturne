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

\page cs_ug_volume_conditions Volume conditions

Volume zone types
=================

Definition of volume regions.

Physical Properties
===================

When the fluid properties are not constant, the user may define the variation laws
in the GUI or in the \ref cs_user_physical_properties user-defined function,
which is called at each time step. In the GUI, in the item **Fluid properties**
under the **Physical properties** section , the variation laws are defined for the
fluid density, viscosity, specific heat, thermal conductivity and scalar dynamic dynamic diffusivity
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
  dynamic diffusivity for the scalars (conductivity <em>λ</em> for the temperature) are considered
  as constant in time and uniform in space, with the \ref cs_fluid_properties_t::cp0
  and the value associated to the field's \ref diffusivity_ref keyword.
  To assign a variable value to <em>C<sub>p</sub></em>, the user **must** specify
  it in the GUI (with a user law) or set \ref cs_fluid_properties_t::icp to
  1 for \ref cs_glob_fluid_properties in \ref cs_user_parameters.c,
  and assign a value for each cell to the array `cpro_cp` which can be
  accessed through \ref CS_F_(cp)->val.

- In the same manner, to have a variable dynamic diffusivity for a given
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

Using physical property libraries
---------------------------------

### CoolProp

When code_saturne has been built with [CoolProp](http://www.coolprop.org)
support, that library may be used to compute fluid properties selected directly
from the GUI.

Note than when properties are variable, CoolProp will default to using
[tabulated data](http://www.coolprop.org/coolprop/Tabular.html) properties
for better performance. The choice of the backend may otherwise
be forced by calling \ref cs_physical_properties_set_coolprop_backend
from the \ref cs_user_model user-defined function.

Note that when using tables, CoolProp stores tabulations in the
`$HOME/.CoolProp/Tables` directory by default, though this can be modified
by defining the `ALTERNATIVE_TABLES_DIRECTORY` environment variable.
Each table requires about 20 Mb in size, and will be created the first time
a given fluid is used in a computation. The initialization time for that
run will reflect this, though subsequent runs will simply read the data.
Tabulation data may be cleaned by removing or cleaning this directory, and its
size may otherwise be limited by using the `MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB`
environment variable.


Initialization
==============

\anchor gui_initialization
\image html gui_initialization.png "Variables initialization"

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

Head losses
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

Source terms
============

momentum, ...
TODO

Basic User source terms {#cs_ug_source_terms}
=======================
*Function called every time step.*

Assume, for example, that the user source terms modify the equation of a
variable \f$\varphi\f$ in the following way:
\f$\rho\frac{\partial \varphi}{\partial t}+\ldots = \ldots + S_{impl}\times\varphi+S_{expl}\f$
The example is valid for a velocity component, for a turbulent variable
(\f$k\f$, \f$\varepsilon$, $R_{ij}\f$, \f$\omega\f$, \f$\varphi\f$ or \f$\overline{f}\f$)
and for a scalar (or for the average of the square of the fluctuations
of a scalar), because the syntax of the \ref cs_user_source_terms function
in the cs_user_source_terms.c file is similar.

For user scalars
----------------

The source terms in the transport equations related to the user scalars
(passive or not, average of the square of the fluctuations of a scalar,
\...) can be filled in thanks to the GUI or the `cs_user_source_terms`
user file. Without the GUI, the `cs_user_source_terms`  function is used
to add source terms to the transport equations related to the user
scalars, this function is called every time step, once for each user
scalar. The user must provide the arrays `st_imp` and `st_exp` related
to each scalar. `st_imp` and `st_exp` must be set to 0 for the scalars
on which it is not wished for the user source term to be applied.

With the GUI
-------------

The GUI can be used if the source terms are proportional to the volume
of the cells or the volume of the fluid in the cell.

By example if we have:

\f$\text{\texttt{st\_exp[c\_id]}}= 0.5*\text{\texttt{cell\_vol[c\_id]}}\f$

\f$\text{\texttt{st\_imp[c\_id]}}= -0.75*\text{\texttt{cell\_vol[c\_id]}}\f$

We can define them directly in the GUI

\anchor gui_user_scal_def_init
\image html gui_user_scal_def_init.png "Transported species and scalars definition"

\anchor gui_select_scalare_source_terms
\image html gui_select_scalare_source_terms.png "select the scalar source terms"

\anchor gui_term_source_meg
\image html gui_source_term_meg.png "initialization of source terms"