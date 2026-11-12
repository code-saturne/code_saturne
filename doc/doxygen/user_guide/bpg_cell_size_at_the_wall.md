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

\page bpg_cell_size_at_the_wall Cell size at the wall

[TOC]

## Hydraulic diameter

The definition of the hydraulic diameter for an internal flow is:
\f[
D_h = \frac{4S}{P}
\f]
with \f$S\f$ standing for the flow surface (cross-section area the fluid flows through)
and \f$P\f$ the wetted perimeter (perimeter of \f$S\f$). For a pipe with a circular
section, the hydraulic diameter is equal to the usual diameter.

## Non-dimensional distance

One generally uses the non-dimensional quantity \f$y^+ = \dfrac{y u^*}{\nu}\f$ to define
the mesh size close to the wall. In this relation:
- \f$y\f$ is the distance from the centre of the mesh cell to the wall (\f$\text{m}\f$);
- \f$u^*\f$ is the friction velocity at the wall (\f$\text{m/s}\f$), representing the wall shear stress \f$\tau_w\f$, defined as:
\f[
u^* = \sqrt{\frac{\tau_w}{\rho}}
\f]
- \f$\nu\f$ is the molecular kinematic viscosity of the fluid (\f$\text{m}^2/\text{s}\f$).

Although the exact calculation of \f$y^+\f$ is not possible *a priori*, one may obtain a reasonable
estimation for internal flows by using a friction velocity of 5% to 10% of the characteristic velocity
value \f$U\f$. The value of 5% corresponds to turbulent flow in a smooth pipe long enough for
turbulence to be fully developed (in this case, Moody's diagram yields approximately \f$u^* \approx \sqrt{0.02 / 8} \cdot U \approx 0.05\,U\f$).

As a conclusion, using a rough approximation, \f$y^+\f$ is of the order of:
\f[
y^+ \approx 0.05\,\frac{yU}{\nu} \quad \text{to} \quad 0.10\,\frac{yU}{\nu}
\f]

## Recommendations

Mesh refinement at the wall should take into account the turbulence model selected.
The following rules can be used as a guide; however, they should be adapted based on the specific configuration.

- **Large Eddy Simulation (LES)**: It is advisable to use a mesh as fine as possible when wall phenomena dominate the flow:
  - *In the direction normal to the wall*: the first cell centre should be such that \f$y^+ < 1\f$, with a very progressive coarsening towards the core of the flow: the ratio of two successive cell sizes should not be larger than 1.1 (i.e. one should use at least 5 to 10 cells below \f$y^+ = 20\f$).
  - *In the main flow direction*: the size \f$\Delta x\f$ of the cells should be such that \f$\Delta x^+ < 40\f$ and the ratio of two successive cell sizes should not be larger than 1.1 or 1.2 (or in general < 2).
  - *In the direction transverse to the flow and parallel to the wall*: the size \f$\Delta z\f$ of the cells should be such that \f$\Delta z^+ < 20\f$ and the ratio of two successive cell sizes should not be larger than 1.1 or 1.2 if no other information is available.

- **v2-f ("phi model") and EBRSM ("low-Reynolds" models)**: It is compulsory to resolve the boundary layer:
  - *In the direction normal to the wall*: it is compulsory that the first cell centre be such that \f$y^+ < 1\f$ and at least 5 to 10 cells should be placed below \f$y^+ = 20\f$.

- **k-omega models (including k-omega SST)**: The cell size at the wall is of little importance:
  - *In the direction normal to the wall*: it is advised, however, to build the mesh so that there are no cells below \f$y^+ = 1\f$ (otherwise the results may become sensitive to the mesh size because of the asymptotic behaviour of the variable \f$\omega\f$).

- **Other RANS models – k-epsilon and RSM – ("high-Reynolds" models)**: Resolution of the viscous sublayer should be avoided:
  - *In the direction normal to the wall*: the mesh should be prepared so that \f$y^+ > 20\f$, preferably around \f$y^+ \approx 30\f$, and it is desirable to build the mesh so that there are at least 5 to 10 cells below \f$y = 0.1\,D_h\f$.

The simulation log file (`run_solver.log` or legacy `listing`) provides useful information to check the adequacy of the mesh for a given model. It is also good practice to use post-processing to visualize the \f$y^+\f$ field that *code_saturne* generates (no attention should be paid to the zero values that may appear locally, as these merely indicate the presence of stagnation, impinging, or detachment regions). These checks are necessary if the user is to adapt or modify the mesh appropriately.

## Origin of the limitations to the wall cell size

The rules proposed above stem from the fact that essentially two approaches may be used to deal with the near-wall region:

- **A fine discretization of the boundary layer**, associated with a physical model capable of representing the viscous sublayer (either a "low-Reynolds" approach or Direct Numerical Simulation (DNS)).
- **A discretization relying on wall functions**, such that the viscous sublayer and buffer layer are entirely contained in the first cell adjacent to the wall by way of an analytical wall law (accounting for the variations of the physical quantities across the boundary layer, most often under the hypothesis that the flow is fully developed and 2D on an infinite flat plate). Thus, the physical model used for the second and subsequent grid cells does not have to model the viscous sublayer.

### Low-Reynolds RANS models

Low-Reynolds RANS models are designed to take into account the physical quantities across the whole boundary layer. However, this requires that the discretization is fine enough to resolve the steep near-wall gradients. This is why one must ensure that \f$y^+\f$ is of the order of 1 at the first cell centre.

### High-Reynolds RANS models

High-Reynolds RANS models, on the contrary, produce incorrect results in the viscous sublayer (\f$y^+ < 10\f$ to \f$20\f$) because they are not designed to resolve boundary layers: wall functions are a solution to this problem, but they are applied only in the cells adjacent to the wall. Hence, it is compulsory to ensure that the rest of the mesh is also outside the viscous sublayer. That is why one tries to ensure that \f$y^+\f$ is of the order of 30 at the centre of the first cell adjacent to the wall.

The use of scalable wall functions for RANS models may alleviate the problems if the local \f$y^+\f$ drops below \f$y^+ < 30\f$. Of course, if the mesh is so fine that several layers of cells are located inside the viscous sublayer (i.e. a "low-Reynolds" mesh), the difficulty remains because the high-Reynolds turbulence equations are not physically valid in the viscous sublayer.

### Mesh topology at the wall (tetrahedra vs. prisms / hexahedra)

A purely tetrahedral grid will cause problems at the wall. In *code_saturne*, boundary conditions are applied at boundary cell faces; however, tetrahedral grids contain cells that have a vertex on the boundary but no face on the boundary. This situation creates local perturbations, as the boundary condition gets applied to cells whose face touches the boundary while not to others whose vertex touches the boundary.

This is one of the main reasons why walls should be meshed with a boundary layer of hexahedra or prisms where it is desirable to capture the boundary layer accurately.

Moreover, it is generally considered that the maximal thickness of a developed boundary layer is approximately \f$0.1\,D_h\f$ for internal flows. If it is important to capture the boundary layer, it should be discretized with several cells: this is why it is advisable to use at least 5 to 10 cells inside the region extending \f$0.1\,D_h\f$ from the wall (from a practical point of view, this is generally a constraint for high-Reynolds models).

### Large Eddy Simulation (LES)

For LES, if the influence of the wall is of primary importance (heat transfer, thermal fatigue, flow-induced vibrations, etc.), one typically wishes to obtain accurate local or transient data at the wall without relying on restrictive wall-function hypotheses. However, resolving the near-wall structures in LES at high Reynolds numbers is often unaffordable due to the prohibitive CPU cost and mesh resolution constraints.

This is why *code_saturne* activates a wall function by default for LES. If the mesh is fine enough (\f$y^+ < 1\f$, \f$\Delta x^+ < 40\f$, \f$\Delta z^+ < 20\f$), the wall function's empirical influence vanishes, and the computation naturally transitions to a wall-resolved LES or quasi-DNS.
