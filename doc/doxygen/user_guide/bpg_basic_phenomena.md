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

\page bpg_basic_phenomena Basic phenomena

## Moving free surface

The movement of the free surface should be considered whenever it is important
to account for:

- transients (for example: pool-filling or draining);
- dispersion of scalar tracers/components by the free surface motion;
- hydrodynamic effects on immersed structures;
- formation of free-surface vortices due to the variation in space or in time of
  the water level (to the knowledge of the authors, no predictive formula exists);
- swell and wave-breaking; let \f$H_h\f$ stands for the crest-to-trough wave height,
  \f$L_h\f$ for the wave length, \f$H\f$ for the water depth: the limit of wave-breaking
  for beaches with a gentle slope is \f$H_h/L_h < 0.14\f$ in deep water
  and \f$H_h/H < 0.8\f$ in shallow water; an approximate generalization of these two
  criteria is \f$H_h/L_h < 0.14 \tanh(2\pi H/L_h)\f$, see (Aelbrecht 2000).

The key non-dimensional numbers to define each case are:

- Froude number \f$F = U/(gH)^{1/2}\f$ with \f$F \gg 1\f$ for supercritical flows and
  \f$F \ll 1\f$ for subcritical flows (the Froude number plays a similar role for
  free surface flows as the Mach number does in aerodynamics)
- Reynolds number \f$Re = UH/\nu\f$ and Weber number \f$We = \rho U^2 H/\sigma\f$
  (if \f$Re > 30\,000\f$ and \f$We > 120\f$, the surface tension has no influence on
  the formation of free-surface vortices)

**In code_saturne:**

The Euler approach with a moving mesh (a.k.a. ALE) is generally well suited for
large deformations of the free surface ("large" with respect to the cell size of
mesh; for example: one may consider that a deformation is large if there are
more than 20 cells for one wave length). Wave-breaking cannot be captured with
this approach.

The displacement of the free-surface may be computed or specified by the user.
If it is computed, the mass conservation may not be perfectly satisfied; this
drawback is associated with the algorithm that is implemented to move the
vertices and the faces of the mesh (and there is no Volume Of Fluid "VOF"
method implemented in code_saturne for the time being).

## Particle laden flows

The particles feature should be used if the following phenomena are to be studied:

- transport, deposition, re-entrainment
- chemistry (wall interaction)

In particular, it is important to decide if:

- the particles follow the fluid or if they should be tracked separately
  because of their own inertia or fluctuating motion (this is especially
  important to capture deposition phenomena)
- the influence of the particles on the flow must be accounted for

A non-dimensional number is important for this configuration:

- The non-dimensional relaxation time-scale of the particles

  \f[
  \tau_p^+ = \frac{(\rho_p/\rho_f) d_p^2 u_*^2}{18 \nu_f^2}
  \f]

  represents the time that is necessary for the velocity of a
  particle to adapt to that of the ambient flow (\f$\rho_p\f$ and \f$\rho_f\f$
  represent the densities of the particles and of the fluid respectively, \f$d_p\f$
  the diameter of the particles, \f$u_*\f$ the wall friction velocity (roughly
  estimated as 5 or 10% of the mean bulk velocity) and \f$\nu_f\f$ the kinematic
  molecular viscosity of the fluid). One may consider that the particles
  follow the fluid if \f$\tau_p^+ \ll 1\f$. Moreover, if \f$\tau_p^+ < 0.1\f$ (generally
  for particles with a diameter smaller that one micron) the deposition is
  affected by Brownian motion. For more details, one may refer to Peirano 2006.

**In code_saturne:**

If it is necessary to account for particles, at least two types of approaches
may be considered:

- Eulerian: represent the particles as a passive scalar tracer, with a
  diffusivity to be defined, under the hypothesis that they follow the fluid
  and that they do not have any influence on the flow itself (that is the choice
  that has been made for the coal combustion module).
- Lagrangian: represent the particles with a statistical approach using a
  Lagrangian method. This method, more general, does not require the hypothesis
  that the particles follow the fluid. If it is necessary, the influence of the
  particles on the flow may be accounted for. The method is useful to study
  deposition and re-entrainment; however, it is important to underline that
  the physical phenomena are complex: for particles whose relaxation
  timescale \f$\tau_p^+\f$ is lower than 1 or 0.1, the standard modelling
  available in code_saturne is not appropriate and more advanced and more
  specific models are required.

## Incompressible, dilatable and compressible flows

By definition, a flow is "compressible" if the fluid density is variable;
if not it is defined as "incompressible"[^1]. For compressible flows,
the pressure waves travel at a finite velocity and the Mach number (\f$M=U/c\f$)
shows the relative importance of the fluid velocity to the compression wave
velocity. In reality, the two types of flow are not clearly distinct since
compressible flows start to behave like incompressible flows when the
Mach number becomes small.

For Mach numbers \f$M=U/c > 0.3\f$, phenomena travelling at the speed of the
sound begin to be important[^2] and it becomes necessary to take into
account the variations of the density due to pressure, temperature, species…
The Boussinesq approximation (that takes into account the density variations
only through the presence of the gravity force in the momentum equation)
is not sufficient any more: the mass equation must contain the unsteady term
(i.e.: \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$) so that
acoustic waves are accounted for.

Only weakly compressible flows at \f$M=U/c < 0.3\f$ will be considered here and
it will be assumed that the effects of the phenomena travelling at the speed
of the sound are negligible (otherwise, code_saturne compressible flow module
should be used). However, the flow may still be compressible, since the density
is not necessarily uniform and constant. Several cases must be considered,
depending on the amplitude of the relative variations of the density.

Before the different cases are presented, the Boussinesq approximation must be
introduced. This linearized approach accounts for the density variations,
when they are small enough, only through the gravity force that appears in the
momentum equation. With this approach, the mass conservation reduces to a steady
constraint on the velocity[^3]: \f$\mathrm{div}(u) = 0\f$.

[^1]: « Dilatable » usually refers to a compressible flow for which the
      variations of the density are assumed to be due to variations of the
      composition or of the temperature (and not to variations of the pressure
      or of the velocity), and for which the possible phenomena associated
      with pressure waves (that are assumed to be infinitely fast) can be
      neglected.

[^2]: The characteristic value \f$M = 0.3\f$ is the usual limit. It comes from
      the following analysis, and ensures that the density variations
      remain "small enough" (see for example (Viollet 1997) or (Wilcox 1997)).
      One considers a fluid accelerating out of a pressurized reservoir:
      the total enthalpy conservation on a streamline (without volume force
      and for steady-state conditions) leads to defining the total enthalpy,
      from which the "total" temperature (or "reservoir temperature") is
      deduced. For a perfect gas, the total temperature reads:
      \f$T_t = T (1+(\gamma-1)/2\, M^2)\f$. Under the hypothesis that the
      entropy does not change
      (\f$P/P_t = (T/T_t)^{\gamma/(\gamma-1)} = (\rho/\rho_t)^\gamma\f$),
      one obtains the "reservoir" pressure
      \f$P_t = P (1+(\gamma-1)/2\, M^2)^{\gamma/(\gamma-1)}\f$ and
      the corresponding density
      \f$\rho_t = \rho (1+(\gamma-1)/2\, M^2)^{1/(\gamma-1)}\f$.
      With these formulae, one can see that the variation of \f$\rho\f$
      remains lower than 2% for \f$M < 0.3\f$.

[^3]: This simplified mass conservation equation comes from dimensional
      analysis considerations, starting from the full mass conservation
      equation \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$,
      with \f$\rho\f$ standing for the density and \f$u\f$ for the fluid velocity.
      The equation can also be written as
      \f$(1/\rho)\, \partial \rho/\partial t + u\, \mathrm{grad}(\rho)/\rho + \mathrm{div}(u) = 0\f$.
      The characteristic length-scale \f$L\f$ and time-scale \f$\delta t\f$
      depends on the physical configuration under consideration. With
      this notation, the magnitude of the three terms are respectively
      \f$(\Delta\rho/\rho)(1/\delta t)\f$, \f$(\Delta\rho/\rho)(u/L)\f$ and \f$u/L\f$.
      For \f$\Delta\rho/\rho \ll 1\f$, the second term is negligible compared to
      the third one; moreover, if \f$\delta t\f$ is not much smaller that \f$L/u\f$
      (the time-scale characteristic of the material waves), the first
      term is also negligible, so that the equation reduces to
      \f$\mathrm{div}(u) = 0\f$.

It is usually considered that the Boussinesq approximation generally applies
if[^4] \f$\Delta\rho/\rho < 0.1\f$ (see for example: (LeQuéré 1992)). Indeed,
(Paillère 2000) presents a configuration with \f$\Delta\rho/\rho = 0.01\f$ for
which the approximation is valid and a configuration with
\f$\Delta\rho/\rho = 0.6\f$ for which the approximation is not valid any more.
Of course, this limit (\f$\Delta\rho/\rho < 0.1\f$) is not absolute:
it merely indicates that the Boussinesq approximation is valid
for "sufficiently small" density variations.

The cases that must be considered are the following, depending on the density
variation magnitude:

- for small variations of the density,
  i.e. \f$\Delta\rho/\rho \ll 1\f$ (usually \f$\Delta\rho/\rho < 0.1\f$),
  the full mass conservation equation
  (\f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$)
  may be used, but one may retain an approximation that is valid under
  some hypotheses:
  - for steady flows (\f$\partial \rho/\partial t=0\f$), the Boussinesq
    approximation may be used;
  - if the flow is not steady, one still may rely on the Boussinesq
    approximation (and use \f$\mathrm{div}(u) = 0\f$ as the mass equation)
    provided the time-scale associated with the variations of the density be
    the convective scale (\f$U/L\f$), which is precisely what is assumed here
    with \f$M < 0.3\f$. Since \f$\Delta\rho/\rho \ll 1\f$, one may also use the
    approximation \f$\mathrm{div}(\rho u) = 0\f$.
- For significant variations of the density (usually \f$\Delta\rho/\rho \geq 0.1\f$),
  the Boussinesq approximation is not valid any more:
  - for steady flows, (\f$\partial \rho/\partial t=0\f$), one may retain
    \f$\mathrm{div}(\rho u) = 0\f$ as the mass conservation equation;
  - if the flow is not steady, it is necessary to use the full mass
    conservation equation \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$.

[^4]: For the air, considered as a perfect gas at atmospheric pressure, the
      limit \f$\Delta\rho/\rho < 0.1\f$ characterizes a temperature variation
      \f$\Delta T/T < 0.1\f$ (with \f$T\f$ in Kelvin), i.e. a variation of 30 K for
      \f$T=300\f$ K.

As a conclusion, for \f$M < 0.3\f$, and under the hypothesis that the effects of
the phenomena travelling at the speed of the sound are negligible
(pressure waves, faster than the material waves), the mass conservation
equation that may be used is as follows:

| \f$M < 0.3\f$ | Steady | Unsteady |
|---|---|---|
| \f$\Delta\rho/\rho < 0.1\f$ | \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$<br>or \f$\mathrm{div}(\rho u) = 0\f$<br>or \f$\mathrm{div}(u) = 0\f$ | \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$<br>or \f$\mathrm{div}(\rho u) = 0\f$<br>or \f$\mathrm{div}(u) = 0\f$ |
| \f$\Delta\rho/\rho \geq 0.1\f$ | \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$<br>or \f$\mathrm{div}(\rho u) = 0\f$ | \f$\partial \rho/\partial t + \mathrm{div}(\rho u) = 0\f$ |

*Table 1: mass conservation equation for M < 0.3*

For more detail, (Gray 1976) provides an example for the derivation of the mass,
momentum and temperature equations where Taylor expansions are used to write the
physical properties (and in particular the density) as functions of temperature
and pressure variations. The hypotheses required for the Boussinesq
approximation to be valid, appear clearly in this systematic approach. For the
specific cases considered in this paper (water and air at 1 atmosphere and 15°C),
the authors obtain a domain of validity defined from the relative variations
of the temperature and of the pressure, of the partial derivatives of the
properties with respect to these two variables, of the Rayleigh and Prandtl
numbers, and of characteristic length- and time-scales.

The non-dimensional numbers that must be evaluated in this configuration are
the following:

- Mach number \f$M = U/c\f$ (\f$M > 0.3\f$: compressible flow)
- relative variation of the density \f$\Delta\rho/\rho\f$
- time-scale characteristic of the flow \f$U/L\f$ and of the boundary conditions

**In code_saturne:**

For \f$M > 0.3\f$ (compressible flows), specific numerical schemes and physical
models are required. A numerical scheme for compressible flows is available
in code_saturne but it has benefited from little feedback and its use requires
some expertise (boundary conditions, thermodynamics for fluids other than
perfect gases, multi-component mixture, …).

For \f$M < 0.3\f$ the standard scheme in code_saturne can be used if
\f$\mathrm{div}(\rho u) = 0\f$ is a valid equation for mass conservation (Table 1),
i.e. except for unsteady flows with \f$\Delta\rho/\rho \geq 0.1\f$ (in this case
(for \f$\Delta\rho/\rho \geq 0.1\f$), a modification of the algorithm of
code_saturne is necessary).

## Natural/forced convection, laminar/turbulent flows

For \f$M < 0.3\f$, the reduced Froude number is
\f$Fr = U/(g\, (\Delta\rho/\rho)\, H)^{1/2}\f$. Natural convection effects
are negligible for \f$Fr \gg 1\f$ (for example for \f$Fr > 10\f$); otherwise, the
gravity force and the density variations must be accounted for (at least using
the Boussinesq approximation, i.e. with the buoyancy term in the momentum
equation).

If gravity has an influence (natural convection), one should evaluate:

- the gradient Richardson number that allows to determine if gravity effects
  inhibit turbulence (\f$Ri > 0.2\f$);
- the Rayleigh number that makes it possible to determine if the regime is
  laminar or turbulent; this is useful to a priori evaluate thermal fluxes
  (using correlations), for comparison to the computational results;

If gravity does not have any influence (forced convection), one should evaluate:

- the Reynolds number that allows us to determine if the regime is laminar
  or turbulent; this is useful to evaluate thermal fluxes a priori (using
  correlations, such as Colburn's) and head losses, for comparison to the
  computational results;

The non-dimensional numbers to evaluate in that case are the following:

- Reduced Froude number \f$Fr = U/(g\, (\Delta\rho/\rho)\, H)^{1/2}\f$
- Gradient Richardson number \f$Ri\f$
- Rayleigh number
- Reynolds number

**In code_saturne:**

When thermal phenomena are neglected or the flow is dominated by forced
convection, it is advised to use:

- High Reynolds number mesh (i.e. a coarse wall grid resolution):
  - the k-epsilon with linear production as the default choice
    (`cs_glob_turb_model->model=CS_TURB_K_EPSILON_LIN_PROD`)
  - the SSG Reynolds Stress Model (`cs_glob_turb_model->model=CS_TURB_RIJ_EPSILON_SSG`),
    whenever secondary motion or turbulent mixing is involved
- Low-Reynolds number mesh (i.e. a fine wall grid resolution):
  - The v2f model as the default choice (`cs_glob_turb_model->model=CS_TURB_V2F_PHI`)

The k-omega SST model (`cs_glob_turb_model->model=CS_TURB_K_OMEGA`) may be
used if it is not clear whether the mesh refinement at the wall makes the
mesh suitable for a high Reynolds or a low Reynolds approach. It is not
advised to use this model otherwise. If possible the v2f model should be used.

For mixed or natural convection phenomena, it is recommended that a
sufficiently fine mesh at the wall should be used (usually possible when
the Rayleigh number is small). If it is not possible to use such a fine
(low Reynolds number) mesh, one may expect poor results for natural or
mixed convection. In this case, there is no really best turbulence model
but it is advised to use k-epsilon or compare the results for a variety
of different turbulence models (eg k-epsilon with linear production
(`cs_glob_turb_model->model=CS_TURB_K_EPSILON_LIN_PROD`), k-omega SST model
(`cs_glob_turb_model->model=CS_TURB_K_OMEGA`) and SSG Reynolds Stress Model
(`cs_glob_turb_model->model=CS_TURB_RIJ_EPSILON_SSG`)).

LES will be used only in cases where local and instantaneous data are required,
or when such phenomena may influence the results that are looked for
(for example: mixing by large structures that RANS models may not be able
to capture). In such cases, the standard Smagorinsky model will be used by
default: for complex configurations, the results produced with this model
are most often as satisfying as those obtained with more advanced models.
The LES model WALE may also be of interest (in code_saturne this model runs
approximately twice as fast as the Smagorinsky model and does not generate
turbulent viscosity in laminar flows which is a well known problem with
the Smagorinsky model).

One may also define the roughness of the wall in the GUI.

For the "scalars" turbulent fluxes (tracers, concentration, temperature,
enthalpy…), one common choice is to
use the high Reynolds number "SGDH" (single gradient diffusion hypothesis)
model. In code_saturne, this model is used with wall-functions and the
Boussinesq approach so that modelling of the turbulent fluxes is proportional
to the gradient of the advected scalar (the coefficient of proportionality is
the ratio of the turbulent viscosity to a turbulent Prandtl or Schmidt number).

This approach suffers limitations, especially for mixed/natural convection or
for low Reynolds number flows (whatever the convection regime). The impact on
the solutions is difficult to quantify a priori. For non-equilibrium anisotropic
flow at high Reynolds numbers, the modelling approaches may encounter
limitations; wall-functions assume that the boundary layer is in equilibrium
and Boussinesq turbulence modelling supposes that the turbulence is isotropic.
More models are available in code_saturne (GGDH – general gradient diffusion
hypothesis -, AFM – algebraic flux model –, DFM – differential flux model -,
to be used in conjunction with the RSM) which take into account the anisotropic
nature of the velocity and the scalar fields.

In all cases, the user must refrain from varying the turbulent Prandtl and
Schmidt numbers with the configuration studied. Such an adjustment would be
equivalent to a case-dependant modification of the diffusion/mixing predicted
by the model: in practice, this would result in user-driven temperature
or concentration results, denying to CFD computations their predictive potential.

## Secondary motion

Identify the elements bound to create secondary motions (corner vortices,
flow structures downstream a bend, vortices in the dead leg of a T-junction,
cyclones separators…) and their potential impact (local exchange coefficient,
perturbation of a stratified flow, source of thermal fluctuations, particle
capture, tracer advection…). These elements should be borne in mind when
choosing the turbulence model (see the previous paragraph).

## Tracers or passive scalars

**In code_saturne:**

### Nature of the scalar

One considers here the "scalars" that represent tracers in code_saturne
(excluding temperature, enthalpy and the other scalar quantities that may
be computed, such as the turbulent variables for example).

code_saturne solves a generic transport equation[^5] for the scalar \f$Y\f$:

\f[
\frac{\partial(\rho Y)}{\partial t} + \mathrm{div}(\rho u Y) = \mathrm{div}(K\, \mathrm{grad}\, Y) + ST_Y
\f]

Hence, the user must verify that this scalar equation is representative of
the problem of interest from a physical point of view.

In the following, one considers the variable \f$s\f$ that stands for the
concentration of a tracer, in kg/m3, (for example, \f$s\f$ may represent
the mass of salt dissolved in the water, per unit of volume of the solution).

The equation that must be solved for \f$s\f$ is the following[^6]
(obtained from considerations on the conservation of the mass of salt):

\f[
\frac{\partial(s)}{\partial t} + \mathrm{div}(us) = DST_{wor}
\f]

However, this is not the equation that code_saturne solves. So \f$s\f$
should not be chosen as a "scalar" in the sense of code_saturne.

A change of variable is required in order to retrieve the equation solved
for "scalars" in code_saturne. If \f$\rho\f$ stands for the density of the
solution (here, salty water), a new variable \f$s'\f$ is defined as:

\f[
s' = s/\rho
\f]

With this definition, \f$s'\f$ is solution of the equation solved by
code_saturne for "scalars":

\f[
\frac{\partial(\rho s')}{\partial t} + \mathrm{div}(\rho u s') = DST_{s'}
\f]

The variable \f$s'\f$, contrary to \f$s\f$, can be selected as a "scalar" in the
sense of code_saturne.

The user may prefer to use a non-dimensional variable instead of \f$s'\f$.
For example, with \f$s'_0\f$ and \f$a_0\f$ standing for constants coefficients,
one may define \f$X\f$ as:

\f[
X = s'/s'_0 - a_0
\f]

With this definition, \f$X\f$ is solution of the equation solved by code_saturne
for the "scalars":

\f[
\frac{\partial(\rho X)}{\partial t} + \mathrm{div}(\rho u X) = DST_X
\f]

The variable \f$X\f$, as \f$s'\f$, can be selected as a "scalar" in the sense
of code_saturne.

[^5]: \f$ST_Y\f$ stands for the source terms that may appear for the scalar \f$Y\f$,
      excluding the molecular diffusion.

[^6]: \f$DST_{wor}\f$ stands for the diffusion and potential source terms for \f$s\f$
      in this equation (for which the density of the fluid solution does not
      explicitly appear on the left-hand-side).

**Remark 1:** it is important to underline that the non-dimensional variable
              standing for the salt concentration and that may be selected
              as a "scalar" in the sense of code_saturne is not \f$s/s_0 - a_0\f$,
              but (with \f$s'_0 = s_0/\rho_0\f$):

\f[
X = \frac{(s/\rho)}{(s_0/\rho_0)} - a_0
\f]

**Remark 2:** the choice of the boundary conditions for the advected scalar
              must be done in accordance with the change of variable that
              defines the selected "scalar" in the sense of code_saturne.
              Usually, this does not lead to any difficulty for Dirichlet
              conditions (for example, the value imposed at the inlet would
              be defined as:
              \f$X_{inlet} = (s_{inlet}/\rho_{inlet}) / (s_0/\rho_0) - a_0\f$).
              However, troublesome situations may arise for non-zero flux
              conditions, if the boundary condition of the density is not
              explicitly defined.

In conclusion, the choice of the "scalars" that can be considered as variables
in code_saturne requires careful attention. For example, if \f$s\f$ is the
concentration of salt in kg/m3 (mass of salt per unit of volume of solution)
and \f$\rho\f$ the density of the salted water in kg/m3, the following variables
can be selected as "scalars" in the sense of code_saturne:

- \f$s' = s/\rho\f$ (mass of salt per unit of volume of solution / density of the
  solution, i.e. mass fraction of the salt)
- \f$C_{ppm} = s' \times 1\,000\,000\f$ (mass fraction of the salt in
  ppm – parts-per-million)
- \f$X = (s/\rho)/(s_0/\rho_0) - a_0\f$ (non-dimensional mass fraction)

### Value of the diffusivity

For a scalar \f$Y\f$, the equation implemented in code_saturne is,
(as indicated above):

\f[
\frac{\partial(\rho Y)}{\partial t} + \mathrm{div}(\rho u Y) = \mathrm{div}(K\, \mathrm{grad}\, Y) + ST_Y
\f]

However, the diffusion term may be expressed as a function of another variable
\f$C = A\,Y\f$.

\f[
\frac{\partial(\rho Y)}{\partial t} + \mathrm{div}(\rho u Y) = \mathrm{div}(K_v\, \mathrm{grad}\, C) + ST_Y
\f]

Under the hypothesis that the variations in space of \f$A\f$ are negligible,
a common approach is to define \f$K\f$ from \f$K_v\f$ as follows (the impact of
this approximation is usually limited, in particular when the turbulent
diffusion, which is modelled, prevails):

\f[
K = K_v A
\f]

For example:

- \f$K = \rho K_v\f$ for \f$C = \rho Y\f$ (with \f$C\f$ a concentration \f$s\f$ in kg/m3 and
  \f$Y\f$ the associated mass fraction \f$s'\f$)
- \f$K = K_v / C_p\f$ for \f$C = Y / C_p\f$ (with \f$C\f$ the temperature and \f$Y\f$ the
  enthalpy for a perfect gas)

A more general approach may be envisaged, under the hypothesis[^7] that \f$C\f$ may
be expressed under the form \f$C(Y)\f$, one gets
\f$\mathrm{grad}\, C = \frac{d(C)}{dY} \mathrm{grad}\, Y\f$, and it is possible
to derive the expression for \f$K\f$ that is necessary to complete the input
data required by code_saturne to solve the equation on \f$Y\f$:

\f[
K = K_v \frac{d(C)}{dY}
\f]

For example:

- In a binary gaseous mixture with \f$M_a\f$ and \f$M_b\f$ the molar masses of the
  components "a" and "b" and with \f$Y\f$ the mass fraction of the component "a",
  the volume fraction of the component "a" reads
  \f$C(Y) = \dfrac{Y/M_a}{Y/M_a+(1-Y)/M_b}\f$ and the diffusivity is
  \f$K = K_v \dfrac{M_a M_b}{(Y M_b+(1-Y)M_a)^2}\f$.
- In a salt water solution for which the density \f$\rho\f$ depends only on
  the mass fraction of salt \f$Y\f$, the mass of salt per unit of volume of
  salted water reads \f$C(Y) = \rho(Y)Y\f$ and the diffusivity of the salt
  is \f$K = K_v (\rho + Y \,d\rho/dY)\f$.

[^7]: \f$C\f$ may not be a function of \f$Y\f$ exclusively: for example, for a
      concentration of salt, when the density also depends on the
      temperature \f$T\f$, \f$C\f$ is a function of \f$Y\f$ and \f$T\f$
      (\f$C = \rho(Y,T)Y = C(Y,T)\f$). In this case, the approach remains
      valid if it is possible to neglect the effects of the transport
      of mass by the diffusion created by the temperature gradient
      (otherwise, some terms depending on the gradient of the temperature
      should be added to the equation of the scalar \f$Y\f$).

## Coupling with the conduction in solids (conjugate gradient with SYRTHES)

When the objective of the study is to determine a thermal load in a
structure, an independent thermal study with boundary conditions representing
the fluid thermal load with correlations may be sufficient. Equally the
computation of the thermal field in the solid is not always compulsory and
predefined correlations may be sufficient to provide boundary conditions for
the fluid calculation. A fully coupled calculation between the fluid and the
solid is required if one is interested in a transient load, local
characteristics, or if the geometry is too complex for reliable correlations
to be available.

The coupling between fluid and solid computations must be taken into account
only if a mutual influence is suspected. This will be the case especially if
the solid is bound to create thermal bridges between fluid zones that
(would otherwise not "see" each other and that) would not mix spontaneously
(for example: the thermally conductive wall of a pipe containing a stably
stratified flow)

- In that case, it is necessary to evaluate the time-scale of the conduction
  in the solid (\f$L_s^2/\lambda_s\f$, ratio of the square of a characteristic
  length of the solid to the thermal conductivity of the solid) and the
  time-scale characteristic of the fluid, i.e. the minimum between at least
  the characteristic convective time-scale \f$U/L\f$ (velocity/characteristic
  length), a characteristic time-scale representative of the possible gravity
  effects over a height \f$H\f$, \f$(H/(g\Delta\rho/\rho))^{1/2}\f$, and a
  characteristic time-scale for turbulence, \f$k/\varepsilon\f$.
- These estimations will show if it is necessary to adopt a specific approach
  to implement a coupled modelling (for example, if the characteristic
  time-scale of the solid if several orders of magnitude larger than that
  of the fluid, an artificial acceleration of the convergence of the solid
  conduction may be employed). These calculations may even help determining
  if the coupling is simply relevant (it may be useless to couple a solid
  that has a very large thermal time-scale to a fluid subjected to very
  rapid variations about a permanent state).

The non-dimensional numbers to evaluate in that case are the following:

Ratio of the time-scales solid/fluid:

\f[
\frac{L_s^2/\lambda_s}{\min\left(U/L,\ (H/(g\Delta\rho/\rho))^{1/2},\ k/\varepsilon\right)}
\f]

## Steady / Unsteady flow

Three types of flows may be considered:

- statistically permanent flows: one will try and model a permanent state
- flows with an inherent unsteadiness: one will compute the time-scales of
  the different phenomena that may play a role. This makes it possible to
  determine whether some of them may be neglected (for example, if they are
  too slow with respect to the others)[^8] or to determine the time necessary
  to reach a possible permanent state. Examples of time-scales:
  - convection: \f$L/U\f$
  - gravity: \f$(H/(g\Delta\rho/\rho))^{1/2}\f$
  - turbulence: \f$k/\varepsilon\f$
  - vortex shedding: \f$L/(US)\f$, where \f$S\f$ stands for the Strouhal number,
    i.e. the non-dimensional frequency of vortex shedding behind an obstacle
    of characteristic size \f$L\f$.
- Flows for which the unsteadiness is driven by the boundary conditions:
  one must compute the time-scales associated with the boundary conditions
  and compare them to the time-scales of the other phenomena. Of course,
  the analysis depends of the boundary conditions considered (explosion,
  valve closing, mass flow rate ramp…).

[^8]: A priori, this analysis is covered by the calculation of non-dimensional
      numbers such a as the Reynolds number, the Froude number, the Mach number,
      that compare the convection, diffusion, gravity and acoustic phenomena.

The (non-) dimensional numbers to evaluate in that case are:

- Fluid time-scales: \f$U/L\f$, \f$(H/(g\Delta\rho/\rho))^{1/2}\f$, \f$k/\varepsilon\f$
- Strouhal number

## Domain of interest

The boundaries of the domain will be defined and their position justified.
This may require using correlations that provide, for example, the minimal
length for a specific type of flow development (pipe flow, jet, mixing layer…)
so as to define the inlet or the outlet locations.

The source and sinks for the momentum (head loss, deviation…) and
temperature (heat) will be defined and the precision with which they shall
be dealt with will be specified.

The inlet and outlet boundary conditions will be defined and the
uncertainties indicated (at least by providing a minimal value,
a maximal value and a probable value). The upstream and downstream regions
of the domain will be described as much as possible (singularities,
head losses of the circuits…).

The surface condition will be indicated (smooth/rough), and the height of
roughness \f$\zeta\f$ will be provided under a non-dimensional form
\f$\zeta^+ = \zeta u_*/\nu\f$ (\f$\zeta^+ < 5\f$: smooth, \f$\zeta^+ > 70\f$: rough).

The non-dimensional numbers to evaluate in that case:

- Non-dimensional size of the surface roughness

## Non-dimensional numbers

| Symbol | Name | Formula | Meaning | Notes |
|---|---|---|---|---|
| \f$F\f$ | Froude | \f$U/(gH)^{1/2}\f$ | inertia / gravity (or velocity / wave speed) | \f$F > 1\f$: supercritical flows, \f$F < 1\f$: subcritical flows. The Froude number plays the same role for free surface flows as the Mach number does for aerodynamics |
| \f$Fr\f$ | Reduced Froude | \f$U/(g(\Delta\rho/\rho)H)^{1/2}\f$ | inertia / reduced gravity | \f$Fr \gg 1\f$: differential gravity forces are negligible with respect to forced convection |
| \f$Gr\f$ | Grashof | \f$g\beta\Delta T L^3/\nu^2\f$ | reduced gravity / viscous effects | Equivalent to the square of a Reynolds number built on a natural convection velocity |
| \f$M\f$ | Mach | \f$U/c\f$ | inertia / wave propagation | \f$M > 1\f$: supersonic, \f$M < 1\f$: subsonic. \f$c^2 = (\partial P/\partial \rho)\vert_s\f$ and for a perfect gas in particular \f$c^2 = \gamma P/\rho\f$. \f$c\f$ is of the order of 300 m/s in air and of 1500 m/s in water (1 bar, 25°C) |
| \f$Nu\f$ | Nusselt | \f$\Phi L/(\lambda \Delta T)\f$ | non-dimensional thermal flux | Correlations such as Colburn, MacAdams (depending on Re, Pr, Ra) |
| \f$Pr\f$ | Prandtl | \f$\nu/a\f$ | viscous effects / conduction effects | 0.6 to 1: gas; 1 to 20: usual liquids; 1000 to 10 000: oils; 0.005 to 0.05: liquid metals |
| \f$Ra\f$ | Rayleigh | \f$g\beta\Delta T L^3/(\nu a)\f$ | \f$Gr\, Pr\f$ | Characteristic of the natural convection regime (laminar below \f$10^5\f$ for example) |
| \f$Re\f$ | Reynolds | \f$UL/\nu\f$ | inertia / viscous effects | \f$Re > 5000\f$: turbulent flow (this limit may be lower for specific types of flow) |
| \f$Ri\f$ | Gradient Richardson | \f$\beta g(\partial T/\partial z)/(2s{:}s)\f$ | gravity / turbulence | \f$Ri > 0.2\f$: turbulence inhibited. pure shear: \f$2s{:}s = \tfrac{1}{2}(\partial U/\partial y)^2\f$; pure impact: \f$2s{:}s = 2(\partial U/\partial x)^2\f$ |
| \f$Sc\f$ | Schmidt | \f$\nu/D\f$ | viscous effects / diffusive effects | Equivalent of the Prandtl number for the diffusivity of the species |
| \f$S\f$ | Strouhal | \f$fL/U\f$ | non-dimensional frequency | - |
| \f$We\f$ | Weber | \f$\rho U^2 L/\sigma\f$ | inertia / surface tension | \f$We \gg 1\f$: surface tension has no influence |

## References

(Aelbrecht 2000) Cours sur la houle, Ecole Nationale des Ponts et Chaussées, D. Aelbrecht, Note EDF R&D HP-72/2000/030/A, 2000

(Gray 1976) The Validity of the Boussinesq Approximation for Liquids and Gases, Donald D. Gray and Aldo Giorgini, International Journal of Heat and Mass Transfer, 19, 545-551, 1976

(LeQuéré 1992) A Chebyshev Collocation Algorithm for 2D Non-Boussinesq Convection, Le Quéré P., Masson R., Perrot P. Journal of Computational Physics 103, 320-335, 1992

(Paillère 2000) Comparison of low Mach number models for natural convection problems, Paillère H., Viozat C., Kumbaro A., Toumi I., Heat and Mass Transfer, 36, 567-573, 2000

(Peirano 2006) Mean-field/PDF numerical approach for polydispersed turbulent two-phase flows, Peirano E., Chibbaro S., Pozorski J. et Minier J.P., Progress in Energy and Combustion Sciences, 32(3): 315-371, 2006

(Viollet 1997) Mécanique des fluides à masse volumique variable, Viollet P.L., Presses des Ponts et Chaussées, 1997

(Wilcox 1997) Basic Fluid Mechanics, Wilcox D.C., DWC Industries, 1997

(Roache 1997) Quantification of Uncertainty in Computational Fluid Dynamics, Annual Rev. Fluid Mech., Vol. 29, pp.123-160, Roache, P.J., 1997
