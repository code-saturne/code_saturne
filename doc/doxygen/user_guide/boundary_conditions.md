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

\page cs_ug_boundary_conditions Boundary conditions

[TOC]

Boundary zone types
===================

Based on the defined zones, the GUI can be used to define boundary
types and conditions:

\anchor gui_bc_definitions
\image html gui_bc_parameters.png "GUI: Base boundary condition definitions"

For advanced cases, user-defined functions are available, as usual.
Many [examples](@ref cs_user_boundary_conditions_examples) are provided.

## Zone-based user-defined function definitions

To define or re-define zone-based boundary condition values,
the \ref cs_user_boundary_conditions_setup (or alternatively
\ref cs_user_finalize_setup) functions may be used.

Note that at walls and when using wall laws (which is the case with most
turbulence models), the boundary values prescribed through Dirichlet
or Robin conditions are not currently directly appplied, but used in conjunction
with the wall model. To force true Dirichlet or exchange conditions,
using the legacy boundary condition settings is still needed.

### Verification of the boundary conditions

The code checks the main compatibilities between the boundary
conditions. In particular, the following rules must be respected:

- If the boundary conditions for the velocity belong to the
  "sliding" type (`icodcl` == 4), the conditions for <em>R<sub>ij</sub>-ε</em>
   must belong to the "symmetry" type (`icodcl`=4), and vice versa.

- If the boundary conditions for the velocity belongs to the "friction" type
  (`icodcl` == 5 or 6), the boundary conditions for the turbulent variables must
  belong to the "friction" type, too.

- If the boundary condition of a scalar belongs to the "friction" type, the
  boundary condition of the velocity must belong to the "friction" type, too.

In case of errors, if the post-processing output is activated (which
is the default setting), a special error output, similar to the mesh
format, is produced in order to better identify and locate boundary condition
definition errors.

Boundary conditions with LES
----------------------------

### Synthetic Eddy Method

The \ref cs_user_les_inflow.cpp user-defined function allows to generate the
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
\ref cs_user_les_inflow.cpp user file.

- \ref cs_user_les_inflow_define (required): define parameters of synthetic turbulence
  at LES inflow.
- \ref cs_user_les_inflow_update (advanced): update of the characteristics of a given synthetic turbulence inlet.
- \ref cs_user_les_inflow_advanced :
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
beyond a threshold value. Given the inlet of size <em>h</em><sup>2</sup> of
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

<!-- ======================================================================= -->

Low-level boundary conditions definitions
=========================================

\subpage low_level_boundary_condition_definitions

<!-- ======================================================================= -->

\page low_level_boundary_condition_definitions Low-level finite-volume boundary-condition definitions

Low-level user-defined functions for boundary conditions
========================================================

For definitions using the legacy finite-volume scheme (i.e. not using CDO),
the \ref cs_user_boundary_conditions (in C) may be used.

For more details about the treatment of boundary conditions, the user
may refer to the [theory](@ref theory) guide and user documentation of the
\ref cs_boundary_conditions_set_coeffs and
\ref cs_boundary_conditions_set_coeffs_turb functions.

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

To access these values, it is recommended to use the `get_val_ext*`,
`get_h_ext*`, and `get_q_ext*` field boundary condition member functions, which
return 1-D or 2-S span-type accessors.

In the case of standard boundary conditions, even without the GUI, it is sufficient
to complete `bc_type[face_id]` and parts of `rcodcl*` arrays, as `icodcl`
and most of `rcodcl*` are filled automatically based on the boundary condition type.
For non-standard boundary conditions, those arrays must be fully completed.

| icodcl  | Description                                                  |
|---------|--------------------------------------------------------------|
| 1       | Prescribed value (Dirichlet condition)                       |
| 2       | Convective boundary condition                                |
| 3       | Flux condition                                               |
| 4       | Symmetry condition                                           |
| 5       | Wall law (friction for velocity, wall value for scalars)     |
| 6       | Rough wall law (deprecated))                                 |
| 11      | Neumann on the normal component, Dirichlet on tangential components. |
| 12      | Imposed value for convection proportional to boundary cell value, imposed flux for diffusion (Robin condition) |
| 13      | Imposed value for convection, imposed flux for diffusion     |
| 14      | Generalized symmetry condition for vectors (Marangoni effect for velocity for example) |
| 15      | User-defined wall law for scalars; for pressure, Neumann for advection and homogeneous Neumann for diffusion (walls with hydro. pressure for compressible module) |

Note that in some cases, a negative value may be used, implying an
automatic conversion:

- For pressure, Dirichlet definition relative to the solved pressure (P*),
  (by default, definition relative to the hydrostatic pressure).
- For enthalpy, conversion from given temperature

Coding of standard boundary conditions {#sec_prg_bc_standard}
--------------------------------------

The standard keywords used by the indicator `bc_type` are, in C:
\ref CS_INLET, \ref CS_SMOOTHWALL, \ref CS_ROUGHWALL, \ref CS_SYMMETRY,
\ref CS_OUTLET, \ref CS_FREE_INLET, \ref CS_FREE_SURFACE,
\ref CS_CONVECTIVE_INLET and \ref CS_UNDEF.

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

- If `bc_type[face_id] == CS_UNDEF`: undefined type face (non-standard case).

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
to `CS_UNDEF` or another value defined by the user (see note at the end
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

  - `f->bc_coeffs->rcodcl1[face_id]` is the value of the variable
    at *I'*, projection of the center of the adjacent cell
    on the straight line perpendicular to the boundary face and crossing its
    center, at the previous time step.

  - `f->bc_coeffs->rcodcl2[face_id]` is CFL number based on the parameter
    <em>C</em>, the distance to the boundary *I'F* and the
    time step: \f$ CFL = \dfrac{C dt }{\centip \centf} \f$.<br><br>

  - `f->bc_coeffs->rcodcl3[face_id]` is not used.

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
    `f->bc_coeffs->rcodcl2[face_id]` must be specified: see [theory](@ref theory).

- If `f->bc_coeffs->icodcl[face_id] == 5`: friction condition, for rough wall
  faces with friction. This condition can not be applied to the pressure.

  - The same rules apply as for a smooth wall, above.

  - The dynamic roughness height is given by
    `CS_F_(vel)->bc_coeffs->rcodcl3[face_id]` only.

  - For the other scalars, the thermal/scalar roughness height is given by
    `f->bc_coeffs->rcodcl3[face_id]` only.

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
for the pressure, an outlet condition (zero Dirichlet if entering, homogeneous Neumann if outgoing) for the velocity, and a
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

Boundary conditions for ALE (Arbitrary Lagrangian Eulerian)  {#sec_ug_bc_ale_legacy}
-----------------------------------------------------------

### Fixed displacement on vertices

For a better precision concerning mesh displacement, one can also assign
values of displacement to certain internal and/or boundary vertices. To do this
one needs to fill the \c depale array and the vertex-based "mesh_displacement"
field values. This property field tracks the total displacement of vertices
relative to their initial position in the mesh.
\c impale[vtx_id] = 1 indicates that the displacement of vertex
\c vtx_id is imposed.

\note Note that the \c impale array is initialized to 0; if its value
      is not modified, corresponding value in \c depale array will be
      ignored.

During the mesh deformation calculation at each time step, the position
of the vertices whose displacement is fixed (i.e. \c impale=1) is not
computed using the value of mesh velocity at the center of corresponding
cell, but directly filled using the values of \c depale.

If the displacement is fixed for all vertices of a boundary face it is not
necessary to prescribe boundary conditions at this face on mesh velocity.
\c bc_coeffs->icodcl and \c bc_coeffs->rcodcl1 values of the "mesh_velocity"
field will be overwritten:
  - \c icodcl is automatically set to 1 (Dirichlet)
  - \c rcodcl1 will be automatically set to face's mean mesh velocity
    value, that is calculated using the "mesh_displacement" field's values.

If a fixed boundary condition (\c ialtyb[face_id]) = CS_BOUNDARY_ALE_FIXED)
is imposed to the face \c face_id, the displacement of each vertex \c vtx_id
belonging to that face is considered
to be fixed, meaning that \c impale[vtx_id] = 1 and
\c displacment[vtx_id][*] = 0 (where mesh_displacement id a pointer to
the "mesh_displacement's" values array).

### Influence on boundary conditions related to fluid velocity

The effect of fluid velocity and ALE modeling on boundary faces that
are declared as walls really depends on the physical nature of this interface.

Indeed when studying an immersed structure the motion of corresponding
boundary faces is the one of the structure, meaning that it leads to
fluid motion. On the other hand when studying a piston the motion of vertices
belonging to lateral boundaries has no physical meaning, therefore it has
no influence on fluid motion.

Whatever the case, the mesh velocity component that is normal to the boundary
face is always taken into account
(\f$ \vect{u}_{fluid} \cdot \vect{n} = \vect{w}_{mesh} \cdot \vect{n} \f$).

The modeling of tangential mesh velocity component differs from one case
to another.

The influence of mesh velocity on boundary conditions for fluid modeling is
managed and modeled in code_saturne as follows (using array names from
\ref cs_user_boundary_conditions_ale):

- If `ale_bc_type[face_id]` = \ref CS_BOUNDARY_ALE_FIXED, the face is motionless;
  mesh velocity equals 0.

- If `ale_bc_type[face_id]` = \ref CS_BOUNDARY_ALE_IMPOSED_VEL, tangential mesh
  velocity is modeled as a sliding wall velocity in fluid boundary conditions
  unless a value for fluid sliding wall velocity has been specified.

- If `ale_bc_type[face_id]` = \ref CS_BOUNDARY_ALE_SLIDING, tangential mesh
  velocity is not taken into account in fluid boundary conditions.

- If `impale[vtx_id] = 1 for all vertices of a boundary face: tangential mesh
  velocity value that has been derived from vertices displacement is modeled as a
  sliding wall velocity in fluid boundary conditions unless a value for fluid
  sliding wall velocity has been specified by USER in code_saturne Interface or
  in 'cs_user_boundary_conditions' function.

Note that mesh velocity has no influence on modeling of
boundary faces with "inlet" or "free outlet" fluid boundary conditions.

For "non standard" conditions, the user has to manage the influence of
ALE boundary conditions (i.e. mesh velocity) on the ones for Navier Stokes
equations(i.e. fluid velocity). (Note that fluid boundary conditions can be
specified in this function.)

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
