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

\page bpg_numerical_parameters Best practice guide: numerical parameters

[TOC]

Introduction
============

The default numerical parameter values of *code_saturne* are chosen so as to provide try to obtain
good precision and resolution of physical phenomena, while maintaining good performance and reasonable robustness,
based on the experience accumulated with the code.

Main numerical parameters
==========================

Gradients computation
---------------------

By default, gradients are computed using a Green-Gauss (finite volume) formulation using a least-squares gradient evaluation of face values with a standard cell neighborhood (see \ref cs_ext_neighborhood_type_t).
For simple reconstructions not requiring a (Finite Volume) conservative gradient a least-squares gradient is used.

The advantages of the Green-Gauss with least-squares face values computation over the iterative Green-Gauss reconstruction are that it has a fixed cost (equivalent to 2 to 3 iterations, whereas the iterative reconstruction typically requires 4 to 5 iterations on a good quality hexahedral mesh, and often more than 10 on tetrahedral meshes (and possibly not converging at all).

Older code_saturne versions (up to 6.0) used the iterative Green-Gauss reconstruction by default.
For some cases, especially those involving anisotropic tensors and good quality meshes,
the iterative algorithm may still have better precision so should be tested in case of issues
with the default.

Also, it is possible to associate different gradient algorithms with different variables,
and even different operation types. The \ref cs_equation_param_t::imrgra values
are those set for each equation using the GUI and are used by default, but
\ref cs_equation_param_t::d_gradient_r is used for diffusion reconstruction, in operators
that do not require finite/volume conservative gradients. By default, a (less costly)
least-squares gradient is used for most of these operators.

In case of robustness issues, switching to a least-squares gradient on an extended
neighborhood often improves the convergence, at the cost of not strictly adhering to
the requirement that some gradients should adhere to a conservative finite volume formulation
(e.g. pressure gradient and velocity gradients used in some turbulence source terms).
In the worst cases, it is necessary to use a full extended neighborhood, though for
some cases, using a less extreme (opposite face venters or optimized heuristic) option
leads to smoother convergence.

Limiters may also be applied to gradient reconstruction, and may be set in the GUI.
Historically, face-based limiters with a clipping factor of 1.5 were used systematically
with least-squares gradients (making it difficult to determine wheher robustness improvements
were due to the gradient reconstruction or the limiter). That setting was decoupled starting
with code_saturne v9.0. In the current version, both gradient limiters (applying a local
multiplier to a gradient used in a reconstruction) and reconstruction limiters (clipping
a reconstructed value to honor the maximum principle using neighbor cell bounds)
are available. Reconstruction limiters should provide a less extreme and better controlled
measure than gradient limiters where applicable.

Using the default settings where possible is recommended, but for some tetrahedral,
polyhedral, or low-quality meshes, using least-squares gradients with limiters
may be necessary. It is recommended to apply such settings only progressively,
variable by variable (at least separately for velocity, pressure, temperature, and
other variables), in the following (weakly recommended) order:
- Apply reconstruction clipping (see flux reconstruction section below).
- Deactivate flux reconstruction at boundaries (see flux reconstruction section below).
- Apply gradient limiter (face-based).
- Switch to a least-squares gradient with extended neighborhood.

In the rare case where a mesh is perfectly orthogonal (such as an axis-aligned Cartesian mesh),
the iterative Green-Gauss algorithm will converge immediately, and be exact in a finite-volume
sense, so is recommended for all variables.

Convective scheme
-----------------

The first-order convective scheme (upwind) is usually more stable and less
accurate than the second-order schemes (centered, SOLU).
- Indeed, the upwind scheme introduces an artificial numerical diffusion that stabilizes the
solution but increases the error. The corresponding physical effect is similar to that of a
diffusion term (second-order derivative in space) proportional to the local velocity and cell size.
- Moreover, the upwind scheme has the advantage of guaranteeing the _maximum principle_:
without source term, the values of a scalar (concentration, temperature) advected by this
scheme remain in the admissible interval, defined by its minimal and maximal initial values –
keeping in mind that this property may be lost because of the discrete treatment of the
diffusion term on arbitrary meshes.

These pure first- and second-order schemes may be combined:
- To increase the stability of the second-order schemes, a “slope test” is activated by default: it
switches locally and without any explicit warning to the upwind scheme wherever oscillations
of the solution are detected. This may reduce the order in space. This slope test may be
switched off (parameter ISSTPC, 0 by default).
- It is possible to select a blended scheme that uses an interpolation between the values that
the pure first and second-order schemes would have produced (formally,
`blencv*second-order+(1-blencv)*first-order`, with `blencv=1.0` by default).

For accurate RANS computations, the following advice may be followed, without claiming
universality:

- For the velocity: a second-order scheme with the default slope test (no maximum principle is
required and it is desirable to give the priority to accuracy).
- For the turbulence: a first-order scheme, usually adequate enough
- For the scalars (concentration, temperature): a second-order scheme with the default slope
test (this choice shall minimize numerical diffusion but does not guarantee the maximum principle.
The solution should be accurate, with little diffusion, but it may exhibit values slightly outside
its physical bounds.

For LES, the same principles apply, but it is recommended to keep the default values.

If the advected variables exhibit oscillations (in regions with high gradient values,
stratification, distorted mesh…), it may be useful to use one of the following solutions, from
the most accurate to the most robust:
- Use the SOLU convection scheme: it is still second-order, but sometimes more stable than the
centered scheme.
- Use a blended convection scheme with a user-defined percentage of upwind (to begin with,
one may prescribe 20% upwind, i.e. `blencv=0.8`).
- Use an upwind scheme for all the advected variables (the result obtained this way will have to
be considered as a first result, still to be improved).

Flux reconstruction
-------------------

The finite volume method relies on the calculation of convective and diffusive fluxes at the
cell faces. These quantities are calculated from values of the variables evaluated at the
orthogonal projections of the cell centers on the straight line normal to the face and
containing its center. By default, the method takes into account the fact that the mesh may not
be orthogonal: the “flux reconstruction” is activated (\ref cs_equation_param_t::ircflu = 1).
If this option is not used,
the simulation is usually more stable but the discretization is not consistent any more
the solution converges as the mesh is refined, but it “sticks” to the mesh and is potentially
wrong: for a diffusion problem, for example, the isolines are parallel to the mesh lines,
independently of the physical phenomenon). Hence, it is generally compulsory to retain the
flux reconstruction option. However, for RANS calculations and especially with
k-epsilon, it has been observed that switching off the flux reconstruction (`ircflu` = 0)
for *k* and *epsilon* could contribute to stabilizing the calculations when difficulties were
encountered on low quality meshes, and this apparently without deteriorating the quality of
the solution for the mean quantities. *A priori*, this also applies to the k-omega and all
the other RANS variants.

Also, as mesh quality issues are often encountered near boundaries, it is possible to deactivate
the flux reconstruction only at boundary faces (\ref cs_equation_param_t::b_diff_flux_rc = 0)
This slightly degrades the spatial convergence order but should not break the scheme's consistency
(as the fraction of the volume where consistency is broken is reduced as the mesh is refined).

Pressure relaxation
-------------------

It is possible to under-relax the pressure increments calculated
within the pressure correction step (by default, no relaxation is used: `relaxv(ipr)=1.0`).
This is useful for distorted meshes and, a priori, it does not affect the quality of
the solution. From previous experience, 30% is a reasonable value for the under-relaxation
parameter (i.e. "Relaxation of pressure increase" = 0.7 instead of the default 1.0 value).
For example, this option may be activated if the calculation fails after several time steps
for which the number of iterations of the pressure iterative solver has been unusually large
or if one suspects that the quality of the mesh is low in a region where large velocity
values develop over a few time steps.

Linear solvers
--------------

### Convergence threshold

The linear solver convergence threshold is already based on an equation-specific
normalization, so the default values of 10^-5 for most variables and 10^-8 for pressure
should usually be safe. Using 10^-5 for pressure is the default for LES, so may be acceptable,
though the correct computation behavior should be verified, comparing results with the
default at least for a number of time steps. This setting may also be used to produce a first
approximation of the solution.

One may also consider this choice if parametric computations are considered. The first
setup may be run with both 10^-5 and 10^-8 threshold for pressure, and it results are
similar enough, it is reasonable to run computations with the other parameters with
a 10^-5 threshold.

### Solver type

By default hybrid, a hybrid Gauss-Seidel solver (or Jacobi on GPU) is used for most
convected variables. With large time steps, the diagonal dominance of the associated
systems may diminsh and convergence degrades, so switching to a Krylov solver
(especially BiCGSTAB, GMRES, or GCR) may be an interesting choice.
Also, is is observed thant in many cases, the default solver converges slowly for the
firs few iterations, then much faster, so that the Krylov solvers are not competitive anymore.
Jacobi and Gauss-Seidel also seem to be positivity-preserving, wich is not the case
for general iterative solvers.

So as to try to use the best solver, a default heuristic switches to GMRES if the solver seems
to converge too slowly (200 iterations with Gauss-Seidel, 400 with Jacobi).
THis is usually quite efficient for the first few iterations, then not needed anymore,
but in cases where the convergence is incomplete (due to ill-conditioned systems),
it is preferrable to force Gauss-Seidel or Jacobi, because when those solvers do not manage to
converge, the intermediate solutions may be closer to the actual solution than those obtained
with un unconverged Krylov solver, and the computation may converge better on later iterations
in some observed cases.

Also, a multigrid solver or preconditioner for convective (or locally convective) systems is
available, but is not activated by default, though it may be very efficient, as its setup cost
is the equivalent of several iterations, so its use is worthwhile only when convergence is
otherwise slow.

Finally, note that when preconditioning an iterative (Krylov-type) solver with a multigrid
solver, only a **flexible** solver should be used, so a Flexible or inexact conjugate gradient
may be used for symmetric systems, or a GCR solver for non-symmetric systems, but non-flexible
variants such as the basic conjugate gradient, GMRES, of BiCGSTAB should be avoided,
as they may exhibit irregular convergence behavior.

Pressure interpolation
----------------------

The pressure gradient is required at the cell centere for the momentum
equation. Hence, because finite volume techniques are used, pressure face-values are necessary.
Those are usually interpolated (centered interpolation) from values at the neighboring cell
centers. This standard approach is appropriate as long as the pressure gradient is reasonably
continuous. When the local variations of the pressure gradient are large, the approach is not
valid any more (a simple centered interpolation cannot account for the fact that the pressure
gradient is significantly different on both sides of the face under consideration). If this
standard approach is used carelessly, the balance of the discrete terms in the momentum
equation is not satisfied anymore and spurious velocities appear (in particular, this
phenomenon may be encountered in the vicinity of a stratification or close to the borders of
a head loss region).
The **Improved pressure interpolation** option (combining `iphydr` and `icalhy` options)
solves the problem through the calculation of an interface pressure taking into account the local
variation of the pressure gradient. Moreover, if the user has defined a head loss region adjacent
to the domain outlet (i.e. in the cells that have an outlet boundary face), `icalhy = 0` should
be set to avoid noticeable perturbations of the pressure boundary condition.

Time-Marching Scheme
====================

Time-marching algorithm and time step value
-------------------------------------------

The “historic algorithm” of code_saturne is a time-marching algorithm, where a time step
*Dt* is specified from the flow characteristics (velocity *U* and kinematic viscosity *ν*
or conductivity *λ*, or their turbulent equivalents) and from the size of the cells *Dx*.

- **To begin with**, an upper limit *Dtmax* for the time step *Dt* shall be evaluated as the
  maximal time-scale value. *Dtmax* is at least as large as:
  * The time that is required for a fluid particle to pass through the whole domain (or the
    time necessary to cross the domain for a signal transmitted by convection or diffusion
    if these phenomena are dominant).
  * The characteristic time associated with potential gravity effects (*H/(g∆ρ/ρ)*) with *H*
    a characteristic height, *∆ρ* the characteristic variation of the density *ρ* and *g*
    the gravity.
- **Ideally**, the time step value *Dt* (`dtref`) is then set on the basis of target Courant
  (*U Dt/Dx*) and Fourier (*ν Dt/Dx*) numbers. One can use target values of 1 and 10
  respectively, in coherence with the default parameters of the code (maximum Courant
  and Fourier numbers). Quite small default values have been selected voluntarily.
- **In practice**, for most of the cases, the value of the time step will be chosen so that the
  Courant number be around 1 to 5 in the main part of the domain of interest. Quite large values
  are acceptable and will not deteriorate the stability of the computation (up to 10 and even 50,
  as long as they are reached only locally). Small values will not endanger the stability of the
  computation either, but the Courant number should not be too low in the important regions (for
  example, not lower than 0.01). Finally, the user should check *a posteriori* that the Fourier
  number computed by the code remains lower than values ranging from 10 to 1000.

### Transient flows

- The time-marching algorithm with a time step constant in time and uniform in space is the
  standard choice in Code_Saturne (`idtvar=0`).
- If the time step value is limited because of a specific part of the transient (for example
  because of the occurrence of a large velocity during a short period of the transient), it is
  possible to let the time step value change in time (`idtvar=1`). The maximum Courant and
  Fourier numbers will be modified to set the desired targets.

### Steady flows

- The time-marching algorithm with a time step uniform in space (`idtvar` = 0 or 1) makes it
  possible to compute a steady solution as the limit of a transient one that has been started from
  an estimated initial state. The convergence may be slow, in particular if the ratio *UDt/Dx*
  varies significantly over the computational domain; in that case, the algorithm with a time step
  which is variable in space and in time may be used (IDTVAR=2). The maximum Courant and
  Fourier numbers will be modified to set the desired targets.

### Warnings

- Convergence
  * With the time-marching algorithm, the calculation must be carried out over
    approximately 5 to 10 times the time required for a fluid particle to cross the whole
    domain (of course, if thermal conduction prevails, the thermal time-scale must also be
    considered).
  * One may also add to the computation the advection of a fictitious passive scalar,
    initially set to 0 over the whole domain and set to 1 at the inlet. The moment when this
    scalar reaches 1 in the whole domain can be used as another convergence indicator
    (not always required, not always sufficient, but generally useful).
  * In all cases, the variation of the quantities provided by the code will be examined (i.e.
    the [time drift](@ref cs_ug_output_time_drift) quantity of the log file associated
    with all the variables but the pressure). One can also use monitoring points to follow
    the evolution of the case-dependent important quantities at user-selected locations
    and at each time step (so that no oscillation may go unnoticed).

- Calculations with heat transfer / scalars
  * It is not advised to use the steady-state algorithm when the flow is influenced by heat
    transfer (through the density variation) or by the distribution of advected scalars
    representing a concentration.
  * On the other hand, if the flow is not seriously affected by the distribution of advected
    scalars, the dynamic features and the scalar transport may be decoupled and it may
    be possible to:
    - Multiply the time step by a factor so as to accelerate the scalar/temperature
      convergence if necessary.
    - Compute the dynamic variables first and deal with the temperature/scalars in
      a second stage on a frozen velocity/pressure field.

- Unsteady RANS
  * The RANS turbulence models are based on statistical decompositions that incorporate the
    turbulent fluctuations into specific variables (k, epsilon, omega, Rij…).
    The use of RANS models for unsteady flows (U-RANS) raises a problem: is there a
    cut-off frequency (and what is its value) that would separate low-frequency
    fluctuations (responsible for variations in time of the mean quantities) from high-frequency
    fluctuations (incorporated into the turbulent variables but potentially unnoticeable on
    the variations in time of the computed quantities)?
  * This is an **open question**. Today, however, unsteady RANS are used without any particular
    precaution. LES are carried out when it is practicable (i.e. when the
    Reynolds number is low enough so that the available computational resources make it possible
    to use a sufficiently fine mesh to resolve enough turbulent structures).

Boundary Conditions
===================

Although the standard boundary conditions generally cover the usual needs, some complements are
provided hereafter.

**Outlet boundary**: for an incompressible flow, the outlet is a theoretical problem per se,
since it is theoretically necessary to know (some characteristics of) the flow downstream of
the outlet to be able to implement the boundary conditions. The standard method (`isolib`)
relies on the hypothesis that the pressure profile does not change in the direction normal to
the outlet. To use this approach in the best possible conditions, it is advised to:
- Select, as much as possible, plane outlets orthogonal to the mean flow so that the hypothesis
  of invariance of the pressure profile in the direction normal to the outlet is as realistic as
  possible.
- Place the outlets sufficiently far away from the regions of interest so as to minimize the
  influence of the boundary conditions on the results of the calculation.
- Place the outlets sufficiently far away from the geometrical perturbations that may produce
  vortices. They may grow arbitrarily when crossing the outlet (at the outlet, if a vortex is
  strong enough, half of it may pump momentum into the domain, while the other half pumps
  momentum out of the domain: even if the flow remains divergence-free, the vortex may
  behave unphysically, sticking to the outlet and growing arbitrarily).
- If it is necessary, move the outlets downstream: the mesh may be extended artificially with
  a few layers of hexahedra or prisms for example.
  * This is *strongly recommended* if the mesh is not otherwise orthogonal at the outlet.
- Set an external Dirichlet value for all the advected scalars (temperature, concentrations):
  it will be used in case the flow re-enters at the outlet (otherwise, a zero-flux condition
  would be used, and make the calculation less stable).
- For large outlets in particular, it may be necessary to set the pressure profile using Dirichlet
  conditions to ensure that the calculation remains stable (for example, a uniform pressure may
  be used: this choice is the responsibility of the user, and depends on how well the pressure is
  known at the outlet).
- For multiple outlets, it is advised to set the pressure at all the outlets but one.
  This approach has the advantage of stabilizing of the calculation, but above all, it is
  physically sound: indeed, it is necessary to provide data characterizing the head loss of
  the circuits that are located downstream of the outlets so as to allow for a physical
  distribution of the flow upstream.

**Wall boundary conditions**: if there is no roughness to specify, the standard wall boundary
conditions do not require any specific attention.

**Inlet conditions**: the standard inlet conditions do not usually require any specific attention.
However, the following advice should be followed:
- Select, if possible, plane inlets orthogonal to the mean flow.
- Place the inlets sufficiently far away from the regions of interest so as to minimize the
  influence of the boundary conditions on the results. It is seldom possible to place the inlet
  sufficiently far away for a fully developed flow to establish (the length that is necessary for
  turbulent a pipe flow to be fully developed is approximately 100 times the hydraulic diameter),
  but it is fortunately seldom useful since the real conditions that should be reproduced in the
  calculation are generally not that of a fully developed flow. However, it is suggested to place
  the inlet approximately at least 10 hydraulic diameters upstream of the regions of interest so
  as to allow for some coherence to develop between the variables that are advected from the
  inlet.
- If it is necessary, displace the inlets upstream: the mesh may be extended artificially with a
  few layers of hexahedra for example.
- For the turbulent variables, provide inlet values that are coherent between each other (the
  default values automatically set by the code are coherent). For example, it would be a
  particularly bad choice to provide a very accurate profile for the turbulent kinetic energy *k*,
  and a value a hundred times too large for the turbulent dissipation *epsilon* (in such a case,
  the data for *k* would immediately be destroyed by the oversized dissipation).
- For LES, it is advised to use the default vortex method or the SEM (Synthetic Eddy Method).
  In practice, for that kind of computation, it is advisable to contact the development team.

Setting-up and checking the computation
=======================================

It is a good idea to check through the topics in this section to avoid common errors that may
lead to erroneous interpretations of the results (the calculation may not fail but instead
produce erroneous results or the user may misinterpret the results).
It is advisable to carry out the verifications suggested here as early as possible
(for example as part of a preliminary computation over the first few steps).

Data input
----------

- Ensure that the mesh has the right units (most often in metres): otherwise, the calculation may
  fail or may seem to progress very slowly – towards wrong results – (for example if a mass flow
  rate designed for an inlet of 1 m is imposed on an inlet of 1 mm or of 1 km ). If a rescaling
  type modification is required, one should keep in mind that Code_Saturne can resize the mesh
  (multiply all dimensions by 1000 for example).
- For inlet conditions, it should be checked right from the start of the calculation that the
  velocity is imposed in the right direction (an error due to the orientation of the axes
  should not be ruled out) and more generally that the variables have the right values and units.
  This may be done through a preliminary post-processing stage, immediately at the beginning of
  the calculation. It should be kept in mind that standard post-processing operates on
  cell-center variables; as a consequence, the exact values set as inlet conditions at the
  boundary faces cannot be observed (however, the difference after a few time steps will
  usually not exceed 1%).
- If free inlet/outlets are supposed to let the fluid out (and not in), it is important to
  ensure right from the start of the calculation that the flow goes in the right direction:
  should it go the wrong way, a flaw may be suspected (for example: a faulty joined mesh,
  an incorrect data input, a time step value too small …)
- When a heat flux is prescribed at the walls, one should double check that the sign of the flux is
  the correct one.
  * In code_saturne, the sign of the flux is positive in the "interior to exterior" direction.
- If variable physical properties have been prescribed (density, viscosity…), the user should plot
  them to check their order of magnitude, their upper and lower limits and the sign of their rate
  of change. Amongst the most common errors, one may encounter the following:
  * The use of a polynomial law outside of its domain of validity (the interval where it
    represents the physical property considered) may lead to large errors, in particular
    with high degree polynomial laws. It is also best to avoid using high degree polynomials,
    replacing them when possible by a series of polynomials of lower degree.
    If *P(T)* represents the density as a function of the
    temperature in the interval *[a ; b]*, it is desirable to compute the numerical values of
    the density as *P(T’)* with *T’ = min(max(T,a),b)*. Ideally, this modification should be a
    mere safety check. Indeed, if *T’* is not exactly equal to *T*, then the user faces a
    problem. On the one hand, if there is no physical reason why the temperature should
    remain in the interval *[a ;b]*, then the choice of the law *P(T)* should be questioned.
    On the other hand, if *T* should physically remain in the interval, then it is the numerical
    scheme that is producing non-physical values; in that case, it is compulsory to use *T’*
    to compute the density, but this modification does not address the substance of the
    problem (which may be associated with the choice of the convective scheme or with
    the quality of the mesh).
  * The coding errors in the user-defined programmes are a common cause of difficulty.
    To avoid these, before starting the computation, one may extract the piece of user
    source code (for example the lines computing *rho(T)*), include it in a separate ad hoc
    programm, let the input data vary artificially (here, the variable *T* would be varied by
    arbitrary steps between its lower and upper limits) and observe the output (here, print
    *rho(T)* and plot it).

Convergence
-----------

Convergence is usually assessed from the evolution of the computational variables at some
well-chosen locations (e.g. monitoring points). For each monitoring point, the code selects the
cell which contains or is closest to the point defined by the user. Hence, the
following advice should be followed:
- Check the position of the cells that have been selected by the code. This can be done
  using the data printed out in the history files (it must not be done using the input data,
  since the latter does not account for the fact that the centers of the selected cells are
  close but not necessarily identical to the user-defined locations). A visual verification
  can be done by positioning each point with Ensight/Paraview on the geometry or by
  plotting with a 2D-plotter the couples (x,y), (y,z) and (x,z) for all the selected cells. For
  example, this preliminary check is particularly important if the quantities of interest of
  the computation is precisely the transient behavior of some variables at some
  monitoring points.
- Record all the variables at all the time steps (or at a high enough frequency) so that no
  oscillation may go unnoticed.
- Check the log file (`run_solver.log`) for the variations of the Courant number: it
  should remain stable (see the values provided previously). It is also advisable to
  visualize the Courant number and its 3D distribution as much as possible.
- Check the log file for the information related to the mass conservation (and to the
  energy conservation, if it is available): the mass conservation should be precisely satisfied,
  with a minimal accuracy of approximately 10^-6 . An inaccurate mass conservation generally
  indicates the imminent failure of the calculation.
- Check the log file for the information related to the energy balance.
  * The **Balance by zone/scalar balance** postprocessing output may be used for this purpose.
  * When the density is constant, the energy conservation should be satisfied with a minimal
    accuracy of approximately 10^-6. When conjugate heat transfer is used with the Syrthes
    code or code_saturne's CDO-based thermal module, the energy balance of code_saturne will
    still be valid with the same accuracy, but the global energy balance (solid+fluid) may not
    be perfect, as coupling may involve interpolation and uses a relaxed explicit time scheme.
    When possible, using code_saturne's internal coupling CHT option will guarantee energy
    conservation.
- Check the log file for the information related to the mean, minimal and maximal
  values of the variables. It is usually a good way to identify the problems of convergence or
  even the input errors.
- Check the log file for the information related to “clippings”: indeed, for the turbulent
  variables (and for the temperature or the scalars, if they are used), an automatic limitation is
  implemented. When the calculation is converged, the number of cells where a variable has
  been clipped should be small (less than 0.1% of the number of cells). If more cells have
  been clipped this may indicate a flaw in the mesh or in the modeling or an unsatisfying level
  of convergence. However, it is possible that the number of cells where a variable has been
  limited be somehow larger (a few %) if a second-order convective scheme is used (especially
  the centered scheme, but also SOLU).
- Check the log file for the “time drift” quantity of the variables: indeed, for a steady-
  state flow, they should diminish by several orders of magnitude between the beginning of the
  calculation and the converged state (by approximately 2 orders of magnitude for the pressure
  and 3 or 4 for the velocity).
- Check the log file for the numbers of iterations of the linear solvers: for a steady-state
  flow, they should diminish quite rapidly and should, for standard computations, tend
  typically towards values of approximately 10 for the velocity and 500 for the pressure (and
  respectively up to 100 and 5000 on large or distorted meshes). These order of magnitude
  values are of course not sufficient to conclude to convergence, but the can help identifying
  quite early a risk of failure of the computation.
