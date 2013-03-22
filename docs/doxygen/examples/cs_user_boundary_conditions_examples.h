/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
*/

/*-----------------------------------------------------------------------------*/

/*!
  \page cs_user_boundary_conditions_examples cs_user_boundary_conditions.f90
 

  \section intro Introduction

  This page provides several examples of code blocks that may be used
  to define boundary conditions in \ref cs_user_boundary_conditions.


  \section cs_user_bc_examples Boundary condition definition examples


  \subsection base_examples Basic examples

  \subsubsection base_loc_var Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet src/user_examples/cs_user_boundary_conditions-advanced.f90 loc_var_dec

  \subsubsection base_init Initialization and finalization

  The following initialization block needs to be added for the following examples:

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 init

  Ad the end of the subroutine, it is recommended to deallocate the work array:

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 finalize

  In theory Fortran 95 deallocates locally-allocated arrays automatically,
  but deallocating arrays in a symmetric manner to their allocation is good
  practice, and it avoids using a different logic for C and Fortran.

  \subsubsection base_body Body

  In the body, we may define several boundary conditions. Here are a few examples.

  \subsubsection base_example_1 Inlet example with hydraulic diameter

  Assign an inlet to boundary faces of group '2' and x < 0.01.

  \warning the <, <=, >, and >= operators may only be used with variables x, y, and z.
  This syntax is not a full equation interpreter, so formulas involving x, y, or z
  are not allowed.

  Set a a Dirichlet value on the three components of \f$ \vect{u} \f$
  on the faces with the selection criterium '2 and x < 0.01' and set a Dirichlet
  to all the scalars \f$ \varia \f$.

  Turbulence example computed using equations valid for a pipe.

  We will be careful to specify a hydraulic diameter adapted
  to the current inlet.

  We will also be careful if necessary to use a more precise
  formula for the dynamic viscosity use in the calculation of
  the Reynolds number (especially if it is variable, it may be
  useful to take the law from \ref usphyv. Here, we use by default
  the 'viscl0" value.

  Regarding the density, we have access to its value at boundary
  faces (romb) so this value is the one used here (specifically,
  it is consistent with the processing in \ref usphyv, in case of
  variable density).

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 example_1

  \subsubsection base_example_2 Inlet example with turbulence intensity

  Assign an inlet to boundary faces of group '3'.

  Set a a Dirichlet value on the three components of \f$ \vect{u} \f$
  on the faces with the selection criterium '3' and set a Dirichlet
  to all the scalars \f$ \varia \f$.

  Turbulence example computed using turbulence intensity data.

  We will be careful to specify a hydraulic diameter adapted
  to the current inlet.

  Calculation of \f$ k \f$ and \f$ \varepsilon \f$
  at the inlet (xkent and xeent) using
  the turbulence intensity and standard laws for a circular pipe
  (their initialization is not needed here but is good practice)

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 example_2

  \subsubsection base_example_3 Assign an outlet to boundary faces of group 'outlet'

  Outlet:
  - zero flux for velocity and temperature, prescribed pressure
  - Note that the pressure will be set to P0 at the first
  - free outlet face (isolib)

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 example_3

  \subsubsection base_example_4 Wall example

  Assign a wall to boundary faces of group '5'.

  Wall:
  - zero flow (zero flux for pressure)
  - friction for velocities (+ turbulent variables)
  - zero flux for scalars

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 example_4

  \subsubsection base_example_5 Rough wall example

  Assign a rough wall to boundary faces of group '7'.

  Wall:
  - zero flow (zero flux for pressure)
  - rough friction for velocities (+ turbulent variables)
  - zero flux for scalars

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 example_5

  \subsubsection base_example_6 Symmetry example

  Assign a symmetry condition to boundary faces of group '4'

  \snippet src/user_examples/cs_user_boundary_conditions-base.f90 example_6


  \subsection channel_inlet Infinite channel inlet

  \subsubsection channel_inlet_loc_var Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet src/user_examples/cs_user_boundary_conditions-auto_inlet_profile.f90 loc_var_dec

  \subsubsection channel_inlet_init Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \subsubsection channel_inlet_example_1 Body

  Here, we define an inlet boundary condition for a very long channel or duct
  with a section matching the boundary faces of group 'INLET'.

  We fix a mean inlet velocity, and use a feedback loop assuming a fixed-point type
  behavior will allow us to reach a state matching that of a very long inlet channel.

  \warning We assume other boundary conditions are defined before this one
  (ideally, using the GUI).
  
  \warning We also assume that the mesh is orthogonal at the inlet, and we are using a RANS
  (not LES) computation.
  to the current inlet.

  For EBRSM of V2f models, to avoid laminarization at the inlet, the initial velocity
  (at the first time step) is divided by 10 on inlet faces adjacent to the boundary, so
  as to ensure a velocity gradient and turbulent production. Otherwise, the initialization
  may fail.

  \snippet src/user_examples/cs_user_boundary_conditions-auto_inlet_profile.f90 example_1


  \subsection advanced_examples Advanced examples

  \subsubsection advanced_loc_var Local variables to be added

  \snippet src/user_examples/cs_user_boundary_conditions-advanced.f90 loc_var_dec

  \subsubsection base_init Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \subsection advanced_ex_1 Example 1

  Example of specific boundary conditions fully defined by the user,
  on the basis of wall conditions selection
  (mass flow computation, specific logging, ...)

  We prescribe for group '1234' a wall, with in addition:
  - a Dirichlet condition on velocity (sliding wall with no-slip condition)
  - a Dirichlet condition on the first scalar.

  \snippet src/user_examples/cs_user_boundary_conditions-advanced.f90 example_1

  \subsection advanced_ex_2 Example 2

  Example of specific boundary conditions fully defined by the user,
  with no definition of a specific type.

  We prescribe at group '5678' a homogeneous Neumann condition for
  all variables.

  \snippet src/user_examples/cs_user_boundary_conditions-advanced.f90 example_2

  \subsection advanced_ex_3 Example 3

  Example of specific boundary conditions fully defined by the user,
  with the definition of a specific type, for example for future
  selection (mass flow computation, specific logging, ...)

  We prescribe for group '6789' a homogeneous Neumann condition for
  all variables, except for the first
  scalar, for which we select a homogeneous Dirichlet.

  \snippet src/user_examples/cs_user_boundary_conditions-advanced.f90 example_3

*/
