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
\ref cs_user_finalize_setup_wrapper) functions may be used.

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

The \ref cs_user_les_inflow.c user-defined function allows to generate the
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
\ref cs_user_les_inflow.c user file.

- \ref cs_user_les_inflow_define (required): define parameters of synthetic turbulence
  at LES inflow.
- \ref cs_user_les_inflow_update (advanced): update of the characteristics of a given synthetic turbulence inlet.
- \ref cs_user_les_inflow_advanced:
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
beyond a threshold value. Given the inlet of size <em>h<up>2</up><em> of
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