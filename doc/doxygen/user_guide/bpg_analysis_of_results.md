\page bpg_analysis_of_results Analysis of the results

It is assumed that the standard verifications have been carried out as recommended in the previous cards (convergence in time and space, coherence of the mesh with the selected turbulence model, and conservation of the quantities that shall be conserved).

In the present document, the physical relevance of the results obtained with the selected modelling will be examined. It is advisable to pay attention to the following points.

## A posteriori verification of the hypotheses

- Verify *a posteriori* that the hypotheses adopted when selecting the models are effectively valid. In particular, evaluate the relevant non-dimensional numbers from the calculation results:
  - Reynolds number;
  - Rayleigh number;
  - Froude number;
  - \f$y^+\f$ for the mesh at the wall;
  - and so on.

## Length and time scales

### Integral length scale

The integral length scale can be evaluated approximately from RANS results as

\f[
L_T = \alpha \frac{k^{3/2}}{\varepsilon},
\f]

where:

- \f$\alpha\f$ ranges approximately from \f$0.1\f$ to \f$0.3\f$;
- \f$k\f$ is the turbulent kinetic energy;
- \f$\varepsilon\f$ is the associated dissipation.

\f$L_T\f$ represents the size of the large structures. In particular, it must be smaller than the characteristic size of the computational domain.

A significantly too large value of \f$L_T\f$ may indicate that the turbulence model is not well suited. If a first-order model has been selected (\f$k\f$-\f$\varepsilon\f$ or \f$k\f$-\f$\omega\f$), one may consider using a second-order model (\f$R_{ij}\f$).

- **Channel flow:** for a hydraulic diameter \f$D_h\f$,

  \f[
L_T = \min(0.42y,\ 0.1D_h),
\f]

  where \f$y\f$ is the distance from the wall.

- **Shear flows (jets or wakes):** as a first approximation [Rodi, 1984], one may consider that \f$L_T\f$ is 10% of the shear-layer width \f$\delta\f$.

  \f$\delta\f$ is defined as the distance between the two points located on both sides of the shear layer where the velocity differs by 1% from the velocity at infinity. For symmetrical flows, \f$\delta\f$ is the distance between the axis of symmetry and the point where the velocity differs by 1% from the velocity at infinity.

### Turbulence time scale

The turbulence time scale can be evaluated approximately from RANS results as

\f[
\frac{k}{\varepsilon}.
\f]

It indicates the lifetime of a turbulent structure. Using \f$k/\varepsilon\f$ and the mean convective velocity \f$U\f$, one may determine the length required for the turbulence to develop:

\f[
U\frac{k}{\varepsilon}.
\f]

## Turbulent variables

### Ratio \f$\nu_t/\nu\f$

For calculations using the high-Reynolds-number \f$k\f$-\f$\varepsilon\f$ model, \f$\nu_t/\nu\f$ is the ratio between the turbulent dynamic viscosity and the molecular dynamic viscosity.

For calculations with wall functions, this ratio should satisfy

\f[
\frac{\nu_t}{\nu} > 10
\f]

at the wall if the first cell adjacent to the wall is large enough, except at singular points such as detachment or reattachment. This condition corresponds approximately to

\f[
y^+ > 25.
\f]

Indeed, in a plane channel flow,

\f[
\nu_t = 0.42yu_*.
\f]

Otherwise, the boundary-layer representation is probably incorrect, and it may be necessary to modify the mesh or at least investigate the relevance of the wall modelling.

## Correlations

### Head losses

When head losses can be evaluated *a priori* from correlations [Idel'cik, 1960], check that the code provides pressure variations in reasonable agreement. A difference of 20% may be considered acceptable.

A significant deviation may indicate that:

- the mesh is not fine enough; or
- the turbulence modelling is unsatisfactory.

The simplified formulae below estimate the pressure loss \f$\Delta P\f$ in several configurations for an incompressible flow of density \f$\rho\f$.

#### Smooth pipe

For a smooth pipe of hydraulic diameter \f$D_h\f$ and length \f$L\f$, with an established flow

\f[
\frac{L}{D_h} > 50 \text{ to } 100,
\f]

and a turbulent flow

\f[
Re = \frac{UD_h}{\nu} > 5000,
\f]

a rough approximation is

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2} = 0.02\frac{L}{D_h}.
\f]

More accurate correlations are

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
= 0.3164Re^{-0.25}\frac{L}{D_h},
\qquad 5000 < Re < 30000,
\f]

and

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
= 0.184Re^{-0.20}\frac{L}{D_h},
\qquad 30000 < Re < 1000000.
\f]

#### Sudden pipe expansion

For a sudden expansion from section \f$S_A\f$, with velocity \f$U_A\f$, to section \f$S_B\f$:

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U_A^2}
= \left(\frac{S_B-S_A}{S_B}\right)^2.
\f]

#### Sudden pipe contraction

For a sudden contraction from section \f$S_A\f$ to section \f$S_B\f$, with velocity \f$U_B\f$:

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U_B^2}
= \frac{1}{2}\left(\frac{S_A-S_B}{S_A}\right).
\f]

#### Smooth pipe bend

For a smooth 90-degree bend with curvature radius \f$R\f$, hydraulic diameter \f$D_h\f$, and a circular or square cross-section, one may remember as a rough approximation that

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
\f]

decreases from \f$1.20\f$ to \f$0.22\f$ as \f$R/D_h\f$ varies from \f$0.5\f$ to \f$1.0\f$, and remains approximately constant as \f$R/D_h\f$ varies from \f$1.0\f$ to \f$5.0\f$. For a 180-degree bend, multiply the pressure loss by approximately \f$1.4\f$.

More precisely [Idel'cik, 1960, p. 192], the pressure loss consists of a bend term, which decreases with \f$R/D_h\f$, and a friction term, which increases with \f$R/D_h\f$.

For a 90-degree bend:

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
= \frac{0.21}{(R/D_h)^{2.5}} + 0.00035\times90\frac{R}{D_h},
\qquad 0.5 < \frac{R}{D_h} < 1.0,
\f]

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
= \frac{0.21}{(R/D_h)^{1/2}} + 0.00035\times90\frac{R}{D_h},
\qquad 1.0 < \frac{R}{D_h}.
\f]

For a 180-degree bend:

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
= 1.4\frac{0.21}{(R/D_h)^{2.5}} + 0.00035\times180\frac{R}{D_h},
\qquad 0.5 < \frac{R}{D_h} < 1.0,
\f]

\f[
\frac{\Delta P}{\tfrac{1}{2}\rho U^2}
= 1.4\frac{0.21}{(R/D_h)^{1/2}} + 0.00035\times180\frac{R}{D_h},
\qquad 1.0 < \frac{R}{D_h}.
\f]

### Heat exchange

When heat exchange can be evaluated from simple correlations [Sacadura, 1980; Taine, 2003], check that the code predicts a reasonable heat flux.

The simplified formulae below estimate the Nusselt number \f$Nu\f$, i.e. the non-dimensional heat flux, for specific configurations of an incompressible flow, where:

- \f$\rho\f$ is the density;
- \f$T_\infty\f$ is the temperature at infinity;
- \f$\lambda\f$ is the thermal conductivity;
- \f$Pr=\nu/a\f$ is the Prandtl number;
- \f$a=\lambda/(\rho C_p)\f$;
- \f$\nu\f$ is the kinematic viscosity;
- \f$\beta=-\frac{1}{\rho}\left.\frac{\partial\rho}{\partial T}\right|_P\f$ is the density variation with temperature at constant pressure;
- \f$g\f$ is the acceleration due to gravity.

#### Forced convection in a smooth pipe

Consider a smooth pipe of length \f$L\f$ and hydraulic diameter \f$D_h\f$, with

\f[
\frac{L}{D_h} > 60,
\f]

for a turbulent flow such that

\f[
Re = \frac{UD_h}{\nu},
\qquad 10000 < Re < 120000,
\f]

and

\f[
0.7 < Pr < 100.
\f]

Fluid properties are evaluated at \f$(T_\infty+T_p)/2\f$. The mean heat flux over length \f$L\f$ is

\f[
\Phi = Nu\,\lambda\frac{T_\infty-T_p}{D_h},
\f]

with the Colburn correlation

\f[
Nu = 0.023Re^{0.8}Pr^{1/3}.
\f]

#### Natural convection on a vertical flat plate

For a vertical flat plate of height \f$L\f$ at constant temperature \f$T_p\f$, with a Rayleigh number

\f[
Ra = \frac{g\beta\Delta T L^3}{\nu a},
\qquad 10^9 < Ra < 10^{13},
\f]

fluid properties are evaluated at \f$(T_\infty+T_p)/2\f$. The mean heat flux is

\f[
\Phi = Nu\,\lambda\frac{T_\infty-T_p}{L},
\f]

with the Mac Adams correlation

\f[
Nu = 0.13Ra^{1/3}.
\f]

### Jets

For an incompressible flow, a round jet emitted from an orifice into a large domain develops in two parts:

1. A **potential-core cone**, in which the velocity remains identical to the velocity at the orifice, progressively spreads over a distance of approximately eight times the orifice diameter.
2. The jet velocity then decreases inversely with the distance from the orifice [Viollet, 1997, citing Abramovitch, 1963].

For the second part, where the centreline velocity decreases, correlations are available [Hug, 1970]. They may help assess the relevance of the results or determine *a priori* the minimum computational-domain size required for the outlet to be sufficiently far away.

Let:

- \f$z\f$ be the distance from the orifice;
- \f$d\f$ be the orifice diameter;
- \f$u_b\f$ be the jet velocity at the orifice;
- \f$c\f$ be the concentration of a tracer, with \f$c=c_b=1\f$ in the jet at the orifice and \f$c=c_a=0\f$ in the ambient fluid.

For constant

\f[
C_p = \left.\frac{\partial h}{\partial T}\right|_P,
\f]

the concentration may represent a non-dimensional temperature:

\f[
c = \frac{T_a-T}{T_a-T_b},
\f]

where \f$T_b\f$ is the temperature in the jet at the orifice and \f$T_a\f$ is the ambient temperature.

The maximum velocity \f$u_m\f$ and maximum concentration \f$c_m\f$ on the jet centreline are

\f[
u_m = 6.2u_b\frac{d}{z},
\f]

\f[
c_m = 5.6\frac{d}{z}.
\f]

### Plane mixing layer

For an incompressible flow, the width \f$e_{cm}\f$ of a mixing layer between two velocities \f$u_1>u_2\f$ is proportional to the distance \f$x\f$ [Viollet, 1997, citing Papamoshkou and Roshko, 1988]:

\f[
e_{cm} = 0.17x\frac{u_1-u_2}{\tfrac{1}{2}(u_1+u_2)}.
\f]

### Backward-facing step

For the flow behind a backward-facing step with a Reynolds number based on the step height of 5100, the DNS of Le and Moin indicates that the reattachment point is located at approximately six times the step height.

See: <http://cfd.mace.manchester.ac.uk/ercoftac/>.

### Vortex shedding behind a circular cylinder

For vortex shedding behind a circular cylinder[^1] of diameter \f$L\f$:

- for an \f$L\f$-based Reynolds number ranging from 200 to 10000, Roshko [1953] measured a Strouhal number of approximately \f$0.21\pm0.01\f$;
- for Reynolds numbers greater than 10000, the Strouhal number increases slightly;
- Roshko [1961] reported measurements of approximately \f$0.27\f$ for a Reynolds number close to \f$5\times10^6\f$.

## References

- **Hug, M.** (1970). *Mécanique des fluides appliquées*. École Nationale des Ponts et Chaussées.
- **Idel'cik, I. E.** (1969). *Mémento des pertes de charges*. Eyrolles, Paris. Cited in the source as “Idel'cik 1960”.
- **Okajima, A.** (1982). “Strouhal Numbers of Rectangular Cylinders”. *Journal of Fluid Mechanics*, **123**, 379-398.
- **Rodi, W.** (1984). *Turbulence Models and Their Application in Hydraulics: A State of the Art Review*. IAHR.
- **Roshko, A.** (1953). *On the Development of Turbulent Wakes from Vortex Streets*. National Advisory Committee for Aeronautics, Technical Note 2913, California Institute of Technology, NACA, Washington.
- **Sacadura, J. F.** (1980). *Initiation aux transferts thermiques*. Tec & Doc.
- **Taine, J., & Petit, J.-P.** (2003). *Transferts thermiques - Introduction aux sciences des transferts* (3rd ed.). Dunod.
- **Viollet, P. L.** (1997). *Mécanique des fluides à masse volumique variable*. Presses des Ponts et Chaussées.

[^1]: For an infinite prismatic obstacle with a square cross-section, the measured values are close to those observed for an infinite circular cylinder. For Reynolds numbers from 100 to 10000, the Strouhal number associated with lift is approximately 0.13 [Okajima, 1982], while the Strouhal number associated with drag - and therefore with vortex shedding - is close to 0.21. The drag frequency is approximately twice the lift frequency because, during one lift period, two vortices detach from the obstacle: one from the upper part and one from the lower part. The drag completes one cycle each time a vortex is shed.
