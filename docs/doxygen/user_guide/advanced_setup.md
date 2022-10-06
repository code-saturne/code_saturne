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

- \subpage low_level_boundary_condition_definitions
- \subpage advanced_specific_physics
- \subpage advanced_coal_and_gas_combution
- \subpage advanced_fuel_oil_combustion
- \subpage advanced_radiative_thermal
- \subpage advanced_conjugate_heat_transfer
- \subpage advanced_particle_tracking
- \subpage advanced_compressible
- \subpage advanced_electric_arcs
- \subpage advanced_coupling
- \subpage advanced_ale
- \subpage advanced_atmospheric
- \subpage advanced_turbomachinery
- \subpage advanced_cavitation

<!-- ----------------------------------------------------------------------- -->

\page low_level_boundary_condition_definitions Low-level finite-volume boundary-condition definitions

Low-level user-defined functions for boundary conditions
========================================================

For definitions using the legacy finite-volume scheme (i.e. not using CDO),
the \ref cs_user_boundary_conditions (in C) may be used. The Fortran equivalent,
\ref cs_f_user_boundary_conditions may still be used but is deprecated, and not
described here (please refer to the documentation from previous versions if needed).

For more details about the treatment of boundary conditions, the user
may refer to the theoretical and computer documentation [@theory] of the
function `condli` (for wall conditions, see `clptur`) (to access this
document on a workstation, use `code_saturne info –guide theory`).

From the user point of view, the boundary conditions are fully defined
by the arrays specifying the boundary type and conditions
assigned to each boundary face.

An array common to all variables defines the boundary condition type relative
to the base flow:

- `bc_type[n_b_faces]`

Other arrays are specific to each solved variable field, and accessible through
its  \ref cs_field_t::bc_coeffs member:

- `icodcl[n_b_faces]` defines the type of boundary condition for the variable.
- `rcodcl1[n_b_faces*dim]`, `rcodcl2[n_b_faces*dim]`, and
  `rcodcl3[n_b_faces*dim]` contain the associated numerical values (value of the
   Dirichlet condition, of the flux \...).

In the case of standard boundary conditions, even without the GUI, it is sufficient
to complete `bc_type[face_id]` and parts of `rcodcl*` arrays, as `icodcl`
and most of `rcodcl*` are filled automatically based on the boundary condition type.
For non-standard boundary conditions, those arrays must be fully completed.

Coding of standard boundary conditions {#sec_prg_bc_standard}
--------------------------------------

The standard keywords used by the indicator `bc_type` are, in C:
\ref CS_INLET, \ref CS_SMOOTHWALL, \ref CS_ROUGHWALL, \ref CS_SYMMETRY,
\ref CS_OUTLET, \ref CS_FREE_INLET, \ref CS_FREE_SURFACE,
\ref CS_CONVECTIVE_INLET and \ref CS_INDEF.

- If `bc_type[face_id] == CS_INLET`: inlet face.

  * Zero-flux condition for pressure and Dirichlet condition for all other
    variables. The value of the Dirichlet condition must be given in
    `f->bc_coeffs->rcodcl1[coo_id*n_b_faces + face_id]` for each solved variable,
    except for the pressure (where `coo_id` is the coordinate id for vector or
    tensor fields, and 0 for scalars).
    The values of `icodcl`, `rcodc2`, and `rcodcl3` and are filled automatically.
    <br><br>

- If `bc_type[face_id] == CS_SMOOTHWALL`: smooth solid wall face, impermeable
  and with friction.

  * the eventual sliding wall velocity of the face is defined by
    `CS_F_(vel)->bc_coeffs->rcodcl1[coo_id*n_b_faces + face_id]`
    (with `coo_id`  0, 1, or 2).
    The initial `rcodcl1` values of  are zero for the three velocity components
    (and therefore need to be specified only if the velocity is not equal to zero).
    \warning The sliding wall velocity must belong to the boundary face
    plane. For safety, the code only uses the projection of this velocity on
    the face. As a consequence, if the velocity specified by the user does
    not belong to the face plane, the wall sliding velocity really taken
    into account will be different.

  * For scalars, two kinds of boundary conditions can be defined:

    - Imposed value at the wall. The user must write:<br>
      `icodcl[face_id]` = 5<br>
      `rcodcl1[face_id]` = imposed value.<br>

    - Imposed flux at the wall. The user must write:<br>
      `icodcl[face_id]` = 3<br>
      `rcodcl3[face_id]` = imposed flux value<br>
      (depending on the variable, the user may refer to the case `icodcl`=3 of
      \ref sec_prg_bc_nonstandard  for the flux definition).<br>

    - If the user does not fill these arrays, the default condition is zero
      flux.<br><br>

- If `bc_type[face_id] == CS_ROUGHWALL`: rough solid wall face, impermeable
  and with friction.

  * The same definitions for scalars and for eventual sliding wall velocity apply
    as for a smooth wall.

  * Roughness is also specified at each boundary face though associated fields:
    `boundary_roughness` for the dynamic roughness<br>
    `boundary_thermal_roughness` for the thermal and scalar roughness.<br>

- If `bc_type[face_id] === CS_SYMMETRY`: symmetry face (or wall without friction).

   * Nothing to specify in `icodcl` and `rcodcl` arrays.<br><br>

- If `bc_type[face_id] == CS_OUTLET`: free outlet face (or more precisely free
  inlet/outlet with forced pressure)

  The pressure is always treated with a Dirichlet condition, calculated with the
  constraint:
  \f$ \frac{\partial }{\partial n}\left(\frac{ \partial P}{\partial \tau}\right) = 0 \f$

  The pressure is set to <em>P<sub>0</sub></em> at the `CS_OUTLET` face closest
  to the reference point defined by `cs_glob_fluid_properties->xyzp0`.
  The pressure calibration is always done on a single face, even if there are
  several outlets.

  If the mass flow is incoming, the velocity is set to zero and a Dirichlet
  condition for the scalars and the turbulent quantities is used (or zero-flux
  condition if no Dirichlet value has been specified).

  If the mass flow is ougoing, zero-flux condition are set for the velocity,
  turbulent quantities, and scalars.

  Nothing is set in `icodcl` or `rcodcl` for the pressure or the velocity.
  An optional Dirichlet condition can be specified for the scalars and turbulent
  quantities.<br><br>

- If `bc_type[face_id] == CS_FREE_INLET`: free outlet or inlet (based on Bernoulli
  relationship) face.

  * If outlet, the equivalent to standard outlet. In case of ingoing flux,
    the Bernoulli relationship which links pressure and velocity is used (see
    the theory guide for more information). An additional head loss modelling
    the outside of the domain can be added by the user.<br><br>

- If `bc_type[face_id] == CS_FREE_SURFACE`: free-surface boundary
   condition.<br><br>

- If `bc_type[face_id] == CS_CONVECTIVE_INLET`: inlet with zero diffusive flux
  for all transported variables (species and velocity).

  * This allows to exactly impose the ingoing flux, without adding a
    diffusive term.<br><br>

- If `bc_type[face_id] == CS_INDEF`: undefined type face (non-standard case).

  * Coding is done in a non-standard way by filling `icodcl`, `rcodcl1`,
    `rcodcl2`, and `rcodcl3` entries associated to that face (see
    \ref sec_prg_bc_nonstandard).<br><br>

\remarks

- Whatever the value of the indicator `bc_type[face_id]`, if the array
  `icodcl[face_id]` is modified for a given variable (*i.e.* filled with a
  non-zero value), the code will not use the default conditions for that
  variable at face `face_id`. It will take into account only the values of
  `icodcl` and `rcodcl` provided by the user (these arrays must then be
  fully completed, like in the non-standard case).

  For instance, for a normal symmetry face where scalar 1 is associated
  with a Dirichlet condition equal to 23.8 (with an infinite exchange
  coefficient) for a given variable field with pointer `f`:

  ```{.C}
  bc_type[face_id] = CS_SYMMETRY;
  f->bc_coeffs->icodcl[face_id]) = 1;
  f->bc_coeffs->rcodcl1[face_id] = 23.8;

  // rcodcl2[face_id] = cs_math_infinite_r is already the default value.
  ```

  The boundary conditions for the other variables are defined automatically.

- The **gradient** boundary conditions in code_saturne boil down to determine
  a value for the current variable <em>Y</em> at the boundary face
  <em>f<sub>b</sub></em>, that is to say \f$ \varia_\fib \f$, value expressed
  as a function of \f$ \varia_{\centip} \f$, value of <em>Y</em> in
  *I'*, projection of the center of the adjacent cell on the straight
  line perpendicular to the boundary face and crossing its center:
  \f[ \varia_\fib=A_{\fib}^g +B_{\fib}^g \varia_{\centip}. \f]

  For a given face, the pair of coefficients \f$ A_{\fib}^g , \, B_{\fib}^g \f$
  may be accessed using the `f->bc_coeffs->a[face_id]` and
  `f->bc_coeffs->b[face_id]` arrays, where the `f` is a pointer to the
  variable's field structure.

  * In the case of a vector or tensor, where `d` represents `f->dim`
    (3 or 6 respectively), `f->bc_coeffs->a[face_id}` is replaced by
    `f->bc_coeffs->a[d*face_id + i]` for coordinate `i` in the expressions
    above, and `f->bc_coeffs->b[face_id}` is replaced by
    `f->bc_coeffs->b[d*d*face_id + d*i + j]`.<br><br>

- The **flux** boundary conditions in code_saturne boil down to determine the
  value of the diffusive flux of the current variable <em>Y</em> at the boundary
  face <em>f<sub>b</sub></em>, that is to say the
  \f$ D_{\ib} \left(K_\fib, \, \varia \right) \f$,
  value expressed as a function of \f$ \varia_{\centip} \f$, value of <em>Y</em>
  in *I'*, projection of the center of the adjacent cell on the
  straight line perpendicular to the boundary face and crossing its center:

  \f[ D_{\ib} \left(K_\fib, \, \varia \right) = A_{\fib}^f +B_{\fib}^f \varia_{\centip}. \f]

  For a given face, the pair of coefficients \f$ A_{\fib}^f , \, B_{\fib}^f \f$
  may be accessed using the `f->bc_coeffs->af[face_id]` and
  `f->bc_coeffs->bf[face_id]` arrays, where the `f` is a pointer to the
  variable's field structure.

  * In the case of a vector or tensor, where `d` represents `f->dim`
    (3 or 6 respectively), `f->bc_coeffs->af[face_id}` is replaced by
    `f->bc_coeffs->af[d*face_id + i]` for coordinate `i` in the expressions
    above, and `f->bc_coeffs->bf[face_id}` is replaced by
    `f->bc_coeffs->bf[d*d*face_id + d*i + j]`.<br><br>

  The **divergence** boundary conditions in code_saturne boil down to
  determining a value for the current variable <em>Y</em> (mainly the Reynolds
  stress components, the divergence \f$ \divv \left(\tens{R} \right) \f$ used in
  the calculation of the momentum equation) at the boundary face
  <em>f<sub>b</sub></em>, that is to say \f$ \varia_\fib \f$, value expressed
  as a function of \f$ \varia_{\centip} \f$, value of <em>Y</em> in
  *I'*, projection of the center of the adjacent cell on the straight
  line perpendicular to the boundary face and crossing its center:
  \f[ \varia_\fib=A_{\fib}^d +B_{\fib}^d \varia_{\centip}. \f]

  For a given face, the pair of coefficients \f$ A_{\fib}^d , \, B_{\fib}^d \f$
  may be accessed using the `f->bc_coeffs->ad[face_id]` and
  `f->bc_coeffs->bd[face_id]` arrays, where the `f` is a pointer to the
  variable's field structure.

  * In the case of a vector or tensor, where `d` represents `f->dim`
    (3 or 6 respectively), `f->bc_coeffs->ad[face_id}` is replaced by
    `f->bc_coeffs->ad[d*face_id + i]` for coordinate `i` in the expressions
    above, and `f->bc_coeffs->bd[face_id}` is replaced by
    `f->bc_coeffs->bd[d*d*face_id + d*i + j]`.<br>

- Caution: to prescribe a flux (nonzero) to Rij, the viscosity to take
  into account is the `molecular_viscosity` field even if the
  `turbulent viscosity` field (ρ.C<sub>μ</sub>/ε).

Coding of non-standard boundary conditions {#sec_prg_bc_nonstandard}
------------------------------------------

If a face does not correspond to a standard type, the user must
completely fill the arrays `bc_type`, `icodcl`, `rcodcl1`, `rcodcl2`,
and `rcodcl3` fr each variable field. `bc_type[face_id]` is then equal
to `CS_INDEF` or another value defined by the user (see note at the end
of \ref sec_prg_bc_standard). The `icodcl` and `rcodcl` arrays must be
filled as follows:

- If `f->bc_coeffs->icodcl[face_id]` == 1: Dirichlet condition.

  * `f->bc_coeffs->rcodcl1[face_id]` is the value of the variable at the
    given face.
    - For vectors and tensors, `f->bc_coeffs->rcodcl1[face_id]` is replaced
      by `f->bc_coeffs->rcodcl1[n_b_faces*i + face_id]` for coordinate `i`.
    - This value has the units of the variable:<br>
      - <em>m/s</em> for the velocity
      - <em>m<sup>2</sup>/s<sup>2</sup></em> for the Reynolds stress
      - <em>m<sup>2</sup>/s<sup>3</sup></em> for the dissipation
      - <em>Pa</em> for the pressure
      - °C for the temperature
      - <em>J.kg <sup>-1</sup></em> for the enthalpy
      - <em>°C<sup>2</sup></em> for temperature fluctuations
      - <em>J<sup>2</sup>.kg<sup>-2</sup></em> for enthalpy fluctuations

  * `f->bc_coeffs->rcodcl2[face_id]` is the value of the exchange coefficient
    between the outside and the fluid for the variable (usually referred to
    as "exterior" exchange coefficient)..
    - An "infinite" value (`rcodcl2[face_id] == cs_math_infinite_r`)
      indicates an ideal transfer between the outside and the fluid (default
      case).
    - It has the following units (defined in such way that when multiplying
      the exchange coefficient by the variable, the given flux has the same
      units as the flux defined below when `icodcl == 3`):
      - <em>kg.m <sup>-2</sup>.s <sup>-1</sup></em> for the velocity
      - <em>kg.m <sup>-2</sup>.s <sup>-1</sup></em> for the Reynolds stress
      - <em>s.m <sup>-1</sup></em> for the pressure
      - <em>W.m <sup>-2</sup>.°C <sup>-1</sup></em> for the temperature
      - <em>kg.m <sup>-2</sup>.s <sup>-1</sup></em> for the enthalpy

  * `f->bc_coeffs->rcodcl3[face_id]` is not used.<br><br>

- If `f->bc_coeffs->icodcl[face_id] == 2`: radiative outlet.

  It reads \f$ \dfrac{\partial \varia }{\partial t} + C \dfrac{\partial \varia}{\partial n} = 0 \f$,
  where <em>C</em> is a to be defined celerity of radiation.

  - `f->bc_coeffs->rcodcl3[face_id]` is not used.

  - `f->bc_coeffs->rcodcl1[face_id]` is the flux value of the variable
    at *I'*, projection of the center of the adjacent cell
    on the straight line perpendicular to the boundary face and crossing its
    center, at the previous time step.

  - `f->bc_coeffs->rcodcl[face_id]` is CFL number based on the parameter
    <em>C</em>, the distance to the boundary *I'F* and the
    time step: \f$ CFL = \dfrac{C dt }{\centip \centf} \f$.<br><br>

- If `f->bc_coeffs->icodcl[face_id] == 3`: flux condition.

  - `f->bc_coeffs->rcodcl1[face_id]` and  `f->bc_coeffs->rcodcl2[face_id]`
    are not used.

  - `f->bc_coeffs->rcodcl3[face_id]` is the flux value of of the variable
    at the wall. This flux is negative if it is a source for the fluid. It
    corresponds to:

    - \f$ -(\lambda_T+C_p\frac{\mu_t}{\sigma_T})\grad T\cdot\vect{n} \f$
      for a temperature (in \f$ W/m^2 \f$)

    - \f$ -(\frac{\lambda_T}{C_p}+\frac{\mu_t}{\sigma_h})\grad h\cdot\vect{n} \f$
      for an enthalpy (in \f$ W/m^2 \f$).

    - \f$ -(\lambda_\varphi+\frac{\mu_t}{\sigma_\varphi})\grad\varphi\cdot\vect{n} \f$
      in the case of another scalar \f$ \varphi \f$ (in \f$ kg.m^{-2}.s^{-1}.[\varphi] \f$,
      where \f$ [\varphi] \f$ are the units of \f$ \varphi \f$).

    - \f$ -\Delta t\ \grad P\cdot\vect{n} \f$ for the pressure
      (in \f$ kg.m^{-2}.s^{-1} \f$).

    - \f$ -(\mu+\mu_t)\grad U_i\cdot\vect{n} \f$ for a velocity component
      (in \f$ kg.m^{-1}.s^{-2} \f$).

    - \f$ -\mu\grad R_{ij}\cdot\vect{n} \f$ for a
      <em>R<sub>ij</sub>-ε</em> tensor
      component (in \f$ W/m^2 \f$).

- If `f->bc_coeffs->icodcl[face_id] == 4`: symmetry condition, for
  symmetry faces or wall faces without friction.

  This condition can only be( used for velocity components
  (\f$ \vect{U}\cdot\vect{n} = 0 \f$) and the <em>R<sub>ij</sub>-ε</em> tensor
  components (for other variables, a zero-flux condition type is usually used).

- If `f->bc_coeffs->icodcl[face_id] == 5`: friction condition, for wall
  faces with friction. This condition can not be applied to the pressure.

  - For the velocity and (if necessary) the turbulent variables, the values
    at the wall are calculated from theoretical profiles. In the case of a
    sliding wall, the three components of the sliding velocity are given by
    `CS_F_(vel)->bc_coeffs->rcodcl1[coo_id*n_b_faces + face_id]`
    (with `coo_id`  0, 1, or 2), the same as for a standard wall.

  - For scalars, the condition `icodcl[face_id] == 5` is similar to
    `icodcl[face_id] == 1`, but with a wall exchange coefficient calculated from
    a theoretical law. Therefore, the values of
    `f->bc_coeffs->rcodcl1[face_id]` and
    `f->bc_coeffs->rcodcl2[face_id]` must be specified: see [@theory].

- If `f->bc_coeffs->icodcl[face_id] == 5`: friction condition, for rough wall
  faces with friction. This condition can not be applied to the pressure.

  - The same rules apply as for a smooth wall, above.

  - The dynamic roughness height is given by
    `CS_F_(vel)->bc_coeffs->rcodcl3[face_id]` only.

  - For the other scalars, the thermal/scalar roughness height is given by
    `f->bc_coeffs->rcodcl3[face_id]` only.

- If `f->bc_coeffs->icodcl[face_id] == 9`: free outlet condition for the
  velocity. This condition is only applicable to velocity components.<br>
  If the mass flow at the face is negative, this condition is equivalent
  to a zero-flux condition.<br>
  If the mass flow at the face is positive, the velocity at the face is
  set to zero (but not the mass flow).<br>
  `rcodcl` is not used.

- If `f->bc_coeffs->icodcl[face_id] == 14`: generalized symmetry boundary
  condition for vectors (Marangoni effect for the velocity for instance).

  This condition is only applicable to vectors and sets a Dirichlet
  boundary condition on the normal component and a Neumann condition on
  the tangential components.<br>
  For each of the 3 components `coo_id`, the required values are:

  `f->bc_coeffs->rcodcl1[coo_id*n_b_faces + face_id]`: Dirichlet value for the
  `coo_id` coordinate.

  `f->bc_coeffs->rcodcl3[coo_id*n_b_faces + face_id]`: flux value for the
  `coo_id` coordinate.

  Therefore, the code automatically computes the boundary condition to
  impose to the normal and to the tangential components.

### Consistency rules summary

In short, following consistency rules between `icodcl` codes for variables with
non-standard boundary conditions:

- Codes for vector or tensor components must be identical
  (the `icodcl` array is handled as a scalar array even for vectors and
  tensors).
- If `icodcl` = 4 for velocity or Rij, it must be 4 for both.
- If `icodcl` = 5 for a scalar or fluctuations, it must be 5 for the velocity.
  * The same rule applies to `icodcl` = 6.
- If `icodcl` = 5 for velocity or turbulence variables, it must be 5 for all
  such variables.
  * The same rule applies to `icodcl` = 6.

Specific cases
--------------

### Outlet faces

A standard `CS_OUTLET` outlet face amounts to a Dirichlet condition (`icodcl == 1`)
for the pressure, a free outlet condition (`icodcl == 9`) for the velocity, and a
Dirichlet condition (`icodcl == 1`) if the user has specified a Dirichlet value
or a zero-flux condition (`icodcl == 3`) for the other variables.

### Boundary condition types for enthalpy

For enthalpy, prescribed values in `CS_F_(h)f->bc_coeffs->rcodcl3` may
be defined using the temperature instead. In this case,
`CS_F_(h)f->bc_coeffs->icodcl` must be replaced by `- CS_F_(h)f->bc_coeffs->icodcl`
(i.e. its sign must bev inverted) to mark the face for automatic
conversion.

### Boundary condition types for compressible flows

For compressible flows, only one of the following boundary condition types
be assigned:

- \ref CS_SMOOTHWALL (standard wall)
- \ref CS_SYMMETRY (standard symmetry)
- \ref CS_ESICF, \ref CS_SSPCF, \ref CS_SOPCF, \ref CS_EPHCF, or \ref CS_EQHCF (inlet/outlet)

Combining approaches
--------------------

Definitions may be based on standard boundary conditions and extended though
non-stand conditions. For example, `bc_type[face_id]` can be set to
`CS_SMOOTHWALL` (whether through the GUI or \ref cs_user_boundary_conditions),
and for a specific variable, the associated `icodcl` and `rcodcl*` arrays
may be modified.

As always, it is recommended to specify only the values which need to be
modified relative to the GUI definitions and default values, so as to
keep user-defined functions concise, readable, and maintainable.

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

\anchor gui_coal_model
\image html gui_coal_model.png "Thermophysical models - Pulverized coal, coal model"

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
or by the user in \ref cs_user_physical_properties for each
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
Initialization of the options of the variables, boundary conditions, initialization of the variables and
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
the \ref user_initialization_compressible for examples of initializations with
the compressible module). Concerning pressure, density, temperature and specific total energy, only 2 variables out of these 4 are independent.
The user may then initialize the desired variable pair
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

\page advanced_electric_arcs Electric arcs module

The electric module is composed of a Joule effect module (\ref CS JOULE EFFECT)
and an electric arcs module (\ref CS ELECTRIC ARCS).

The Joule effect module is designed to take into account that effect (for instance in glass
furnaces) with real or complex potential in the enthalpy equation. The Laplace forces are not
taken into account in the impulse momentum equation. Specific boundary conditions can be
applied to account for the coupled effect of transformers (offset) in glass furnaces.

The electric arcs module is designed to take into account the Joule effect (only with real
potential) in the enthalpy equation. The Laplace forces are taken into account in the
impulse momentum equation.

Activating the electric arcs module
===================================

The electric arcs module is activated either:

 - in the Graphical User Interface _GUI_: __Calculation features__ --> __Electrical arcs__, the user can choose between _Joule Effect_ for _joule model_ and _Joule Effect_ and _Laplace Forces_ for electric arc
 - or in the user function \ref cs_user_model in cs_user_parameters.c file, by setting the \ref cs_glob_physical_model_flag[\ref CS_ELECTRIC_ARCS] or \ref cs_glob_physical_model_flag[\ref CS_JOULE_EFFECT] parameter to a non-null value.

Initialization of the variables
===============================

The function \re cs_user_initialization allows the user to initialize some of the specific physics variables prompted via \ref cs_user_model. It is called only during the initialization of the calculation. As usual,the user has access to many geometric variables so that the zones can be treated separately if needed (see [Electric arcs example](@ref user_initialization_electric_arcs)).

The values of potential and its constituents are initialized if required.

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
 and can be used for the initialization of the variables in
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

\warning
 In the case of alternating current, attention should be paid to the values of potential
 imposed at the limits: the variable named "real potential" represents an affective
 value if the current is in single phase, and a "real part" if not.

 - For the Joule studies, a complex potential is sometimes needed
 (*in the GUI Electrical model -> three-phase*): this is the  case in particular where the current
 has three phases. To have access to the phase of the potential, and not just to its
 amplitude.

 - For the Joule studies in which one does not have access to the phases, the real
 potential (imaginary part =0) will suffice (*in the GUI Electrical model -> AC/DC*): this is
 obviously the case with
 continuous current, but also with single phase alternative current. In code_saturne
 there is only 1 variable for the potential,  called "real potential". Pay attention to
 the fact that in alternate current, the "real potential" represents a effective value
 of potential, \f$\frac{1}{\sqrt{2}}\,Pp_{max}\f$ (in continuous current there is no
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

\page advanced_coupling Coupling with other domains

code saturne-code saturne coupling
==================================

The user function \ref cs_user_saturne_coupling in \ref cs_user_coupling.c is
used to couple *code_saturne* with itself.

Such couplings allow explicit exchange of boundary conditions and source terms,
possibly with different model choices in each domain.

In the case of a coupling between two *code_saturne* instances, first argument *saturne_name*
of the function \ref cs_sat_coupling_define is ignored.
In case of multiple couplings, a coupling will be matched with available *code_saturne*
instances based on that argument, which should match the directory name for the
given coupled domain see [examples](@ref cs_user_coupling_h_cs_user_saturne_coupling).

Fluid-Structure external coupling
=================================

The function \ref usaste belongs to the module dedicated to external
Fluid-Structure coupling with *code_aster*. Here one defines the boundary
faces coupled with *code_aster* and the fluid forces components which are
given to structural calculation. When using external coupling with *code_aster*,
structure numbers necessarily need to be negative;\n
the references of coupled faces being (*i.e. -1, -2*), etc.
For examples on the function we can see [examples](@ref cs_user_fluid_structure_interaction_h_usaste)

<!-- ----------------------------------------------------------------------- -->

\page advanced_ale ALE (Arbitrary Lagrangian Eulerian) module

Setting options
===============

The **ALE module** may be activated through the Graphical User Interface (GUI)
in the **Calculation features** section. It can also be activated
in the \ref cs_user_model function in \ref cs_user_parameters.c.
See [ALE activation](@ref cs_user_parameters_h_cs_user_ale) for examples.


When activated in the GUI, a **Deformable mesh** page appears in the GUI,
Providing additional options. The user must choose the
type of mesh viscosity and describe its spatial distribution,
see [Ini-ale](@ref gui_ale_mei).

\anchor gui_ale_mei
\image html gui_ale_mei.png "Thermophysical models - mobile mesh (ALE method)"

The options may also be set and completed with the \ref usstr1 function
found in \ref cs_user_fluid_structure_interaction.f90. It is possible to
specify informations for the structure module, such as the index of the
structure (\ref idfstr), and the initial displacement, velocity and acceleration
values (\ref xstr0, \ref xstreq and \ref vstr0). For examples,
see [examples](@ref cs_user_fluid_structure_interaction_h_usstr1).

Mesh velocity boundary conditions
=================================

These boundary conditions can be managed through the GUI, or using the
\ref cs_user_boundary_conditions_ale  function in
\ref cs_user_boundary_conditions.c file.

With the GUI, when the item **Deformable mesh** is selected
in **Calculation features**, specific boundary condition types
for the mobile mesh may be defined for each zone.
The user can choose a boundary condition type for ALE (internal coupling),
including displacements based on a mass-spring model, such as
[CL-ale1](@ref gui_ale_internal_bc).

\anchor gui_ale_internal_bc
\image html gui_ale_internal_bc.png "Mesh boundary conditions - internal mass-spring model coupling"

When at least one boundary is of the "internal coupling" (mass-spring model)
type, a **Coupling parameters** entry appears under
the one for **Boundary conditions**.

\anchor gui_ale_internal_param
\image html gui_ale_internal_param.png "Coupling parameters - internal coupling"

Instead of or in addition to settings with the GUI, the
\ref cs_user_boundary_conditions_ale fonction may be used with the ALE module.
It is used in a similar way to the \ref cs_user_boundary_conditions
in the framework of standard calculations, that is to say array values
are defined for each face to specify the detailed mesh boundary conditions.
See [Examples of boundary conditions for ALE](@ref example_ale2).

Modification of the mesh viscosity
==================================

With the ALE module, the \ref cs_user_physical_properties user-defined function
allows modifying the mesh viscosity.
It is first called before the time loop, and before reading restart files
(so the mesh is always in its initial position at this stage).
The user can modify mesh viscosity values to prevent cells and nodes from huge
displacements in awkward areas, such as boundary layer for example.

Note that for more complex settings, the mesh viscosity could be modified in
\ref cs_user_initialization or \ref cs_user_extra_operations.
The matching field's name is **mesh_viscosity**.

Fluid - structure internal coupling
===================================

In the file \ref cs_user_fluid_structure_interaction.f90 the user provides the parameters of two functions.

- \ref usstr1 is called at the beginning of the calculation. It is used to
  define and initialize the internal structures where fluid-Structure coupling
  occurs. For each boundary face, **idfstr** is the index of the structure the
  face belongs to (if idfstr(ifac) = 0, the face ifac doesn't belong to any
  structure). When using internal coupling, the structure index must be strictly
  positive and smaller than the number of structures.
  The number of *internal* structures is automatically defined with the maximum
  value of the **idfstr** table, meaning that internal structure numbers must
  be defined sequentially with positive values, beginning with integer value
  **1**
  [see here for more details](@ref cs_user_fluid_structure_interaction_h_usstr1).

- The second function, \ref usstr2, is called at each iteration. It is used to
  define structural parameters (considered as potentially time dependent):
  mass m \ref xmstru, friction coefficients  \ref xcstru, and stiffness
  k \ref xkstru.
  The \ref forstr array defines fluid stresses acting on each internal
  structure
  [see here for more details](@ref cs_user_fluid_structure_interaction_h_usstr2).
  Moreover it is also possible to take external forces (gravity for example)
  into account.

<!-- ----------------------------------------------------------------------- -->

\page advanced_atmospheric Atmospheric flows module

[TOC]

Data files
==========

When using the atmospheric module, a file called **meteo** may be added to
a case's **DATA** directorys in order to provide vertical
profiles of the main variables.

Atmospheric domain mesh requirements
====================================

An atmospheric mesh has the following specific features:

 - The boundary located at the top of the domain should be a plane.
So, horizontal wind speed at a given altitude can be prescribed at the top
face as an inlet boundary.

 - Cells may have very different sizes, from very small (near ground or
buildings) to very large (near the top of domain or far from zone of interest).

 - Vertical resolution: from tiny cells (e.g. \f$\Delta\f$\upshape z = 1 m) near
the ground to a few hundreds of meters at the top.

 - Horizontal resolution: from a few meters to hundreds of meters.

 - The length ratio between two adjacent cells (in each direction) should
preferably be between \f$0.7\f$ and \f$1.3\f$.

 - The z axis represents the vertical axis.

A topography map can be used to generate a mesh. In this case, the preprocessor
 mode is particularly useful to check the quality of the mesh (run type Mesh
quality criteria).

Atmospheric flow model and steady/unsteady algorithm
====================================================

The GUI may be used to enable the atmospheric flow module and set up the
following calculation parameters in the **Thermophysical models - Calculation features**
page see [fig:steady](@ref gui_atmospheric_user_s_guide_v92).

The atmospheric flow model
--------------------------

The user can choose one of the following atmospheric flow models:
 - **Constant density**: To simulate neutral atmosphere.
 - **Dry atmosphere**: To simulate dry, thermally-stratified atmospheric flows (enables *Potential temperature* as thermal model).
 - **Humid atmosphere**: To simulate thermally stratified atmospheric flows (air-water mixture) with phase changes
   (enables *Liquid potential temperature* as thermal model). The model is described in \cite Bouzereau:2004 and \cite Bouzereau:2007.

Allowed time-stepping options
-----------------------------

- The pseudo-steady time-stepping is usually chosen. It sets a time step
  variable in space and time. It can be selected if constant boundary
  conditions are used, and usually provides fastest and smoothest convergence.
- The unsteady time-stepping algorithm must be used in presence of time varying
  boundary conditions or source terms (the time step can then be variable in time or constant).

\anchor gui_atmospheric_user_s_guide_v92
\image html gui_atmospheric_user_s_guide_v92.png "Selection of atmospheric model"

\anchor gui_atmospheric_user_s_guide_v93
\image html gui_atmospheric_user_s_guide_v93.png "Selection of steady/unsteady flow algorithm"

\warning
The following points have to be considered when setting the parameters
described above:

  - The potential temperature thermal model and the liquid potential
temperature one (see the paragraph **Atmospheric main variables** for the
definition) requires that the vertical component of the gravity is set to
\f$g_z=-9.81 m.s^{-2}\f$ (\f$g_x=g_y=0 m.s^{-2}\f$),
otherwise pressure and density won't be correctly computed.
  - As well, the use of scalar with drift for atmospheric dispersion requires
the gravity to be set to \f$g_z=-9.81\f$ (\f$g_x=g_y=0 m.s^{-2}\f$), even if the density
is constant.

Physical properties
===================

The specific heat value has to be set to the atmospheric value
\f$C_{p}=1005 J/kg/K\f$.

| **Parameters** | **Constant** \n **density** | **Dry atmosphere** | **Humid atmosphere** | **Explanation** |
|----------------|-------------------------|--------------------|----------------------|-----------------|
|pressure boundary \n condition|Neumann first \n order|Extrapolation|Extrapolation|In case of **Extrapolaion**, \n the pressure gradient is assumed (and set) constant, whereas in case of\n **Neumann first order**, the pressure gradient is assumed (and set) to zero.|
|Improved pressure|no|yes|yes|If yes, exact balance between the hydrostatic part of the pressure gradient and the gravity term \f$\rho\f$g is numerically ensured.|
|Gravity (gravity is assumed aligned with the z-axis)|\f$g_z=0\f$ or \f$g_z=-9.81 m.s^{-2}\f$ (the latter is useful for scalar with drift)|\f$g_z=-9.81 m.s^{-2}\f$|\f$g_z=-9.81 m.s^{-2}\f$|              |
|Thermal variable|no|potential temperature|liquid potential temperature|                        |
|Others variables|no|no|total water content, droplets number|                                   |

Boundary and initial conditions
===============================

The *meteo* file can be used to define initial conditions for the
different fields and to set up the inlet boundary conditions. For the velocity
field, **code_saturne** can automatically detect if the boundary is an inlet boundary or an
outflow boundary, according to the wind speed components given in the
*meteo* file with respect to the boundary face orientation. This is often
used for the lateral boundaries of the atmospheric domain, especially if the
profile is evolving in time. In the case of inlet flow, the data given in the
*meteo* file will be used as the input data (Dirichlet boundary condition)
for velocity, temperature, humidity and turbulent variables. In the case of
outflow, a Neumann boundary condition is automatically imposed (except for the
pressure). The unit of temperature in the *meteo* file is the degree
Celsius whereas the unit in the GUI is the kelvin.

To be taken into account, the *meteo* file has to be selected in the GUI
(*Atmospheric flows* page, see [fig:meteo](@ref gui_atmo_read)) and the check
box on the side ticked. This file gives the profiles of prognostic atmospheric
variables containing one or a list of time stamps. The file has to be put in the
**DATA** directory.

\anchor gui_atmo_read
\image html gui_atmo_read.png "Selection of the *meteo* file"

An example of file *meteo* is given in the directory
data/user/meteo. The file format has to be strictly respected.
The horizontal coordinates are not used at the present time (except when
boundary conditions are based on several meteorological vertical profiles)
and the vertical profiles are defined with the altitude above sea level. The
highest altitude of the profile should be above the top of the simulation domain
and the lowest altitude of the profile should be below or equal to the lowest
level of the simulation domain. The line at the end of the *meteo* file
should not be empty.

If the boundary conditions are variable in time, the vertical profiles for
the different time stamps have to be written sequentially in the *meteo*
file.

You can also set the profiles of atmospheric variables directly in the GUI.
The following boundary conditions can be selected in the GUI:
 - Inlet/Outlet is automatically calculated for lateral boundaries (e.g. North, West\textellipsis ) of the computational domain
(see [fig:inlet](@ref gui_atmospheric_user_s_guide_v95)).
 - Inlet for the top of the domain (see [fig:top](@ref gui_atmospheric_user_s_guide_v96)).
 - Rough wall for building walls (see [fig:walls](@ref gui_atmospheric_user_s_guide_v97)) or for
the ground (see [fig:ground](@ref gui_atmospheric_user_s_guide_v98)).
The user has to enter the roughness length. In case of variable roughness
length, the user has to provide the land use data and the association
between the roughness length values and land use categories.

\anchor gui_atmospheric_user_s_guide_v95
\image html gui_atmospheric_user_s_guide_v95.png "Selection of automatic inlet/ outlet for boundary conditions"

\anchor gui_atmospheric_user_s_guide_v96
\image html gui_atmospheric_user_s_guide_v96.png "Selection of the boundary condition for the top of the domain"

\anchor gui_atmospheric_user_s_guide_v97
\image html gui_atmospheric_user_s_guide_v97.png "Selection of the boundary condition for building walls"

\anchor gui_atmospheric_user_s_guide_v98
\image html gui_atmospheric_user_s_guide_v98.png "Selection of the boundary condition for the ground"

\remark
    If a meteorological file is given, it is used by default to
initialize the variables. If a meteorological file is not given, the user can
use the standard **code_saturne** initial and boundary conditions set up but has to be aware
that even small inconsistencies can create very large buoyancy forces and
spurious circulations.

Boundary conditions based on several meteorological vertical profiles
---------------------------------------------------------------------

In some cases, especially when outputs of a mesoscale model are used, you
need to build input boundary conditions from several meteorological vertical
wind profiles. Cressman interpolation is then used to create the boundary
conditions. The following files need to be put in the \texttt{DATA} directory:
\item All *meteo* files giving the different vertical profiles of
prognostic variables (wind, temperature, turbulent kinetic energy and
dissipation).
 - A file called imbrication_files_list.txt which is a list
of the *meteo* files used.
 -  A separate *meteo* file which is used for the initial conditions
and to impose inlet boundary conditions for the variables for which Cressman
interpolation is not used (for example: temperature, turbulent kinetic energy).
This file must follow the rules indicated previously.

The following files should be put in the SRC directory:
  - The user source file cs_user_parameters.f90. In this file, set
the cressman_flag of each variable, for which the Cressman
interpolation should be enabled, to *true*.

User-defined functions
======================

User-defined functions may be used when the graphical user interface is not
sufficient to set up the calculation. We provide some examples of user file for
atmospheric application:
 - cs_user_source_terms.c: to add a source term in the prognostic equations for forest
   canopy modelling, wind turbine wake modelling... [examples](@ref user_source_terms)
 - cs_user_parameters.f90: to activate the Cressman interpolation.
   For example, it is used to impose inhomogeneous boundary conditions.
   [examples](@ref cs_f_user_parameters_h_usati1)
 - cs_user_extra_operations.c to generate vertical profiles for post processing.
   [examples](@ref cs_user_extra_operations_examples_mean_profiles)
 - cs_user_boundary_conditions.f90: showq how to set up the boundary conditions and to set
   a heterogeneous roughness length... [examples](@ref atmospheric_examples)

Physical models
===============

Atmospheric dispersion of pollutants
------------------------------------

To simulate the atmospheric dispersion of pollutant, one first need to define
the source(s) term(s). That is to say the location i.e. the list of cells or
boundary faces, the total air flow, the emitted mass fraction of pollutant,
the emission temperature and the speed with the associated turbulent parameters.
The mass fraction of pollutant is simulated through a user added scalar that
could be a *scalar with drift* if wanted (aerosols for example).

The simulations can be done using 2 different methods:
 - Prescribing a boundary condition code **total imposed mass flux** for
some boundary faces using the cs_user_boundary_conditions.f90 user function.
 - Using a scalar source term. In this case, the air inflow is not taken
   into account. The user has to add an explicit part to the equations
   for the scalar through the cs_user_source_terms.c file. This is
   done by selecting the cells and adding the source term \ref st_exp
   which equals to the air flux multiplied by the mass fraction, while the
   implicit part \ref st_imp is set to zero.

With the first method, the same problem of sources interactions appears, and
moreover standard Dirichlet conditions should not be used (use
itypfb=i_convective_inlet and icodcl=13 instead) as
the exact emission rate cannot be prescribed because the diffusive part
(usually negligible) cannot be quantified. Additionally, it requires that
the boundary faces of the emission are explicitly represented in the mesh.\n

Finally the second method does not take into account the jet effect of the
emission and so must be used only if it is sure that the emission does not
modify the flow.\n

Whatever solution is chosen, the mass conservation should be verified by using
for example the [cs_user_extra_operations-scalar_balance.c](@ref cs_user_extra_operations_examples_scalar_balance_p) file.

Soil/atmosphere interaction model
---------------------------------

This model is based on the force restore model (\cite Deardorff:1978).
It takes into account heat and humidity exchanges between the ground and the
atmosphere at daily scale and the time evolution of ground surface temperature
and humidity. Surface temperature is calculated with a prognostic equation
whereas a 2-layers model is used to compute surface humidity.

The parameter \ref iatsoil in the file \ref atini0.f90 needs to be equal to one to
activate the model. Then, the source file \ref solvar.f90 is used.

Three variables need to be initialized in the file \ref atini0.f90: deep soil
temperature, surface temperature and humidity.

The user needs to give the values of the model constants in the file
\ref solcat.f90: roughness length, albedo, emissivity...

In case of a 3D simulation domain, land use data has to be provided for the domain.
Values of model constants for the land use categories have also to be
provided.

Radiative model (1D)
--------------------

The 1D-radiative model calculates the radiative exchange between different
atmospheric layers and the surface radiative fluxes.

The radiative exchange is computed separately for two wave lengths intervals

- Calculation in the infrared spectral domain (file \ref rayir.f90)
- Calculation in the spectral range of solar radiation (file \ref rayso.f90)

This 1D-radiative model is needed if the soil/atmosphere interaction model
is activated.

This model is activated if the parameter \ref iatra1 is equal to one in the
file cs_users_parameters.f90.

Atmospheric main variables
==========================

For more details on the topic of atmospheric boundary layers, see \cite stull:1988.

- Definition of the potential temperature:
\f[
\theta =T\left(\frac{P}{P_{r}}\right)^{-\frac{R_{d}}{C_{p}}}
\f]
- Definition of liquid potential temperature:
\f[
\theta_{l} = \theta \left( 1-\frac{L}{C_{p}T} q_{l} \right)
\f]
- Definition of virtual temperature:
\f[
T_{v} = \left(1+0.61q\right)T
\f]
- Gas law:
\f[
P = \rho \frac{R}{M_{d}}\left(1+0,61q\right)T
\f]
with \f$R=R_{d} M_{d}\f$.
- Hydrostatic state:
\f[
\frac{\partial P}{\partial z} = -\rho g
\f]

| **Constant** **name** | **Symbol** | **Values** | **Unit** |
|-----------------|----------|----------|--------|
|Gravity acceleration at sea level|\f$g\f$| \f$9.81\f$|\f$m.s^{-2}\f$|
|Effective Molecular Mass for dry air| \f$M_{d}\f$ | \f$28.97\f$| \f$kg.kmol^{-1}\f$|
|Standard reference pressure | \f$P_{r}\f$ | \f$10^{5}\f$ | \f$Pa\f$|
|Universal gas constant | \f$R\f$| \f$8.3143\f$ | \f$J.K^{-1}.mol\f$|
|Gas constant for dry air| \f$R_{d}\f$ | \f$287\f$ | \f$J.kg^{-1}.K^{-1}\f$|


| **Variable** **name** | **Symbol** |
|---------------------|----------|
|Specific heat capacity of dry air | \f$C_{p}\f$|
|Atmospheric pressure | \f$P\f$|
|Specific humidity | \f$q\f$|
|Specific content for liquid water | \f$q_{l}\f$|
|Temperature | \f$T\f$|
|Virtual temperature | \f$T_{v}\f$|
|Potential temperature | \f$\theta\f$|
|Liquid potential temperature | \f$\theta_{l}\f$|
|Latent heat of vaporization|\f$L\f$|
|Density | \f$\rho \f$|
|Altitude | \f$z\f$|

Recommendations
===============

This part is a list of recommendations for atmospheric numerical simulations.

- Enough probes at different vertical levels in the domain should be used
  to check the convergence of the calculation.
- An inflow boundary condition at the top level of the domain should be set
  (symmetry and automatic inlet/outlet are not appropriate).
- A Courant number too small or too big has to be avoided (see code_saturne Best
  Practice Guidelines). That is the reason why the
  **variable time step in space and in time** option is recommended for steady
  simulations when there are large differences of cell size inside the domain
  (which is generally the case for atmospheric simulations). With this option,
  it can be necessary to change the reference time step and the time step maximal
  increase (by default, the time step increase rate is \f$10\f$).

In some cases, results can be improved with the following modifications:

- In some case, the turbulent eddy viscosity can drop to unrealistically low
  values (especially with \f$k-\varepsilon\f$ model in stable atmospheric condition).
  In those cases, it is suggested to put an artificial molecular viscosity around
  \f$0.1 m^{2}.s^{-1}\f$.
- If the main direction of wind is parallel to the boundary of your computing
  domain, try to set symmetry boundary conditions for the lateral boundaries to
  avoid inflow and outflow on the same boundary zone (side of your domain).
  Another possibility is to use a cylindrical mesh.
- To avoid inflow and outflow on the same boundary zone (side of your domain),
  avoid the case of vertical profile in the input data \texttt{meteo} file with
  changes of the sign of velocity of wind (\f$V_x\f$ or/and \f$V_y\f$).

<!-- ----------------------------------------------------------------------- -->

\page advanced_turbomachinery Turbomachinery module

Introduction
============

Two classical models are available in **code_saturne** for rotor/stator
interactions modelling in turbomachinery computations: the steady
approach which is based on the so-called *Frozen Rotor* modelling
and the *transient rotor/stator* approach which is based on a
sliding mesh technique.

\warning

This section describes these functionalities based on
a single **code_saturne** computation. An alternative rotor/stator coupling based
on coupling of boundary conditions is also possible (and only briefly
described in this section) but it is not recommended.

Meshing recommendations
=======================

Periodicity
-----------

The rotational periodicity treatment is possible only in *Frozen
Rotor*. However, the interface plane between rotor and stator
must match in the azimutal $\theta$ direction:
  - \f$\theta_{min}^{rotor}(z)=\theta_{min}^{stator}(z),\quad\theta_{max}^{rotor}(z)=\theta_{max}^{stator}(z)\f$

for all \f$z\f$ through the rotation axis direction.

Rotor/stator interface
----------------------

Unsteady rotor/stator: in the input mesh(es), the
interface between rotor and stator domains has to be composed of
\underline boundary faces. Then the interface boundary faces are joined
during the computation and become internal faces, as is usual for
mesh joining in the preprocessing stage. A simple way to ensure
joining is not done prematurely is to provide
\underline separated meshes for each rotor or stator domain.
  - *Frozen Rotor*: the interface can be composed of boundary
   faces (in which case the interface boundary faces are joined at
   the beginning of the computation) or of internal faces.

Meshing of the interface region
-------------------------------

As mentioned above, when a rotor/stator interface boundary exists (in
particular for the *unsteady rotor/stator* model), boundary faces
are joined by the solver during the computation, based on the current
rotor position. It is thus important to be aware that the success of
a joining operation is strongly dependant on the
\underline quality of the mesh at the interface. More precisely,
the refinement must be as similar as possible at both sides of the
interface. Moreover, it is reminded that the tolerance parameter of
a joining is a fraction of the shortest edge linked with a vertex of
a joined face. Consequently, cells with high aspect ratios where the
refinement in the azimutal \f$\theta\f$ direction is much coarser than
those in one of the two others can also lead to a joining failure.
In particular, the user should be careful to avoid elongated
viscous layer type cells in curved areas such as a rotor-stator interface.

If the meshes at both sides of the interface are very different
such that the joining fails, advanced joining parameters are
available. However, modifying the mesh is more likely to
succeed. The introduction of a somekind of buffer cells layer on
both sides of the interface should be very valuable. Ideally, each
of the two layers should have the same refinement and a constant
azimutal step (this latter recommandation is relevant only for
*unsteady rotor/stator* model).

Alternative rotor/stator coupling
---------------------------------

If the meshes at both sides of the interface are very different and
can not be modified, a fallback solution is to use the rotor/stator model
based on the boundary conditions coupling.

\Warning: Contrarily to the mesh joining approach, the
boundary conditions coupling approach is not fully conservative.

Turbomachinery dedicated postprocessing functions
=================================================

Useful postprocessing functions relative to the machinery
characteristics are available: postprocessing of the couple on the
rotor walls and postprocessing of the head generated by the machinery.

Data setting, keywords and examples
===================================

Data setting, keywords and examples for turbomachinery computations
(mesh joining or boundary conditions coupling), are provided in [the dedicated doxygen documentation.](@ref turbomachinery)

<!-- ----------------------------------------------------------------------- -->

\page advanced_cavitation Cavitation module

The cavitation module is based on an homogeneous mixture model. The
physical properties (density and dynamic viscosity) of the mixture
depends on a resolved void fraction and constant reference properties
of the liquid phase and the gas phase.

For a description of the user management of the cavitation module,
please refer to [the dedicated doxygen documentation.](@ref cavit)
