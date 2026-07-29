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

\page cs_ug_base_setup Base modeling setup

The documentation for this section is in the process of migration from the pdf
documentation. It is also recommended to check the
[pdf user's guide](../../../user/user.pdf) for sections which may not have been migrated yet.

[TOC]

- \subpage base_setup_initialization
- \subpage base_setup_mass_injection
- \subpage gui_user_law_editor
- \subpage base_var_modify_end_time_step
- [CDO schemes](@ref cs_ug_cdo_hho_base)

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
