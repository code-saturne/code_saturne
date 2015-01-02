/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
  \page parameters Input of calculation parameters (C functions)

  \section cs_user_parameters_h_intro Introduction

  C user functions for definition of model options and calculation parameters.
    These subroutines are called in all cases.

  If the Code_Saturne GUI is used, this file is not required (but may be
    used to override parameters entered through the GUI, and to set
    parameters not accessible through the GUI).

  Several functions are present in the file, each destined to defined
    specific parameters.


  \section cs_user_parameters_h_cs_user_model  Base model related options

  Definition of user variables or properties should be defined here,
  if not already done throught the GUI.

  \section cs_user_parameters_h_cs_user_linear_solvers  Linear solver related options

  By default, Code_Saturne will use a multigrid algorithm for pressure
  and iterative solver for other variables. For a given case, checking
  the setup file resulting from a first calculation will provide
  more info.

  Available solvers include a variety of iterative linear solvers,
  described in more detail at \ref cs_sles_it_create, and a multigrid
  solver, whose definition and settings are described at
  \ref cs_multigrid_create, \ref cs_multigrid_set_coarsening_options,
  and \ref cs_multigrid_set_solver_options.

  Simple options may be set using the GUI, but for more advanced settings
  are described in this section. It is also recommended to read
  the documentation of \ref cs_sles.c (which is a solver definition
  "container"), \ref cs_sles_it.c (iterative solvers, with available
  types \ref cs_sles_it_type_t), and \ref cs_multigrid.c
  (which are actual solver implementations). The API provided
  is extensible, so it is possible for a user to define other solvers
  or link to external solver libraries using this system,
  without requiring any modification to non-user source files.

  The examples which follow illustrate mostly simple setting changes
  which may be useful.

  \subsection cs_user_parameters_h_sles_ex_1 Example: distance to wall

  By default, the wall distance (active only with turbulence models which
  require it) is computed with a preconditionned conjugate gradient.
  The following example shows how to use a multigrid solver for this
  quantity (useful especially if computed repeatedly, such as for ALE).

  \snippet cs_user_parameters-linear_solvers.c sles_wall_dist

  \subsection cs_user_parameters_h_sles_user_1 Example: user variable

  The following example shows how to set the linear solver for a given
  user variable field so as to use a BiCGStab solver with polynomial
  preconditioning of degree 1.

  \snippet cs_user_parameters-linear_solvers.c sles_user_1

  \subsection cs_user_parameters_h_sles_user_1 Changing the verbosity

  By default, a linear solver uses the same verbosity as its matching variable,
  and is not verbose for non-variable quantities. The verbosity
  may be specifically set for linear system resolution, as shown in
  the following example:

  \snippet cs_user_parameters-linear_solvers.c sles_verbosity_1

  \subsection cs_user_parameters_h_sles_mgp_1 Example: advanced multigrid settings

  The following example shows how to set advanced settings for the
  multigrid solver used for the pressure solution.

  \snippet cs_user_parameters-linear_solvers.c sles_mgp_1

  \subsection cs_user_parameters_h_sles_mg_parall Multigrid parallel settings

  In parallel, grids may optionally be merged across neigboring ranks
  when their local size becomes too small. This tends to deepen
  the grid hierarchy, as some parallel rank boundaries are removed.
  Depending on the architecture and processor/network performance
  ratios, this may increase or decrease performance.

  \snippet cs_user_parameters-linear_solvers.c sles_mg_parall

  \subsection cs_user_parameters_h_sles_rad_dom Example: DOM radiation settings

  For DOM radiation models, 1 solver is assigned for each direction
  this allows using a specific ordering for each direction for the
  default Block Gauss-Seidel solver.

  The example below shows how to set a non-default linear solver for
  DOM radiation. Here, we assume a quadrature with 32 directions
  is used (if more solvers than directions are specified, the extra
  definitions will be unused, but this causes no further issues).

  \snippet cs_user_parameters-linear_solvers.c sles_rad_dom_1

  \section cs_user_parameters_h_cs_user_moments  Time moment related options

  Code_Saturne allows the calculation of temporal means or variances,
  either of expressions evaluated through a user function, or
  of expressions of the type \f$<f_1*f_2...*f_n>\f$. The variables
  may be fields or field components. This is done calling either
  through the GUI, or in the user function \ref cs_user_time_moment.
  Each temporal mean is declared using either
  \ref cs_time_moment_define_by_func, or \ref cs_time_moment_define_by_field_ids.

  For each time moment, a starting time step or value is defined. If the
  starting time step number is negative, the time value is used instead.

  The moment values are stored as fields, and updated at each time step,
  using recursive formulas. Before the matching moment computation
  starting time step, a moment's values are uniformly 0.
  For visualization an interpretation reasons, only fields of dimension
  1, 3, 6, or 9 (scalars, vectors, or tensors of rank 2) are allowed, so
  moment definitions not matching this constraint should be split.

  To count defined moments, use the \cs_time_moments_n_moments function,
  whether from Fortran or C. To access the matching fields, use
  \ref time_moment_field_id in Fortran, or \ref cs_time_moment_get_field
  in C.

  \section cs_user_parameters_h_examples Examples

  \subsection cs_user_parameters_h_example_1 Example 1

  In the following example, we define a moment for the mean velocity.
  All components are used (component -1 means all components),
  so the moment is a vector.

  \snippet cs_user_parameters-time_moments.c tmom_u

  \subsection cs_user_parameters_h_example_2 Example 2

  In the next example, we multiply the expression by the density.
  As the density is of dimension 1, and the velocity of dimension 3,
  the resulting moment is of dimension 3.

  \snippet cs_user_parameters-time_moments.c tmom_rho_u

  \subsection cs_user_parameters_h_example_3 Example 3

  In the next example, we define a product of several field components,
  all of dimension 1, as we consider only the x and y components of the
  velocity; for the density, we cas use either component 0 or -1 (all),
  since the field is scalar.

  This moment's computation is also restarted at each time step.

  \snippet cs_user_parameters-time_moments.c tmom_rho_u_v

  \subsection cs_user_parameters_h_example_4 Example 4

  This next example illustrates the use of user-defined functions
  to evaluate expressions. Here, we compute the moment of the sum
  ot two variables (which obviously cannot be expressed as a product),
  so we first need to define an appropriate function, matching the
  signature of a \ref cs_time_moment_data_t function.
  We can name that function as we choose, so naming for clarity is recommmended.
  Note that in this case, the input argument is not used. This argument
  may be useful to pass data to the function, or distinguish between calls
  to a same function.

  Note also that we compute both means and variances here.

  \snippet cs_user_parameters-time_moments.c tmom_simple_sum_data

  In \ref cs_user_time_moments, we can not assign that function to a moments
  definition:

  \snippet cs_user_parameters-time_moments.c tmom_simple_sum

  \subsection cs_user_parameters_h_example_5 Example 5

  In this last example, we compute components of the mean velocity
  in the case of a rotating mesh. As the mesh orientation changes
  at each time step, it is necessary to compensate for this
  rotation when computing the mean, relative to a given mesh position.
  When using the matching moment, it will also be necessary to
  account for the mesh position.

  Here, the same function will be called for each component, so
  an input array is defined, with a different key (here a simple
  character) used for each call.

  \snippet cs_user_parameters-time_moments.c tmom_velocity_rotation_data

  Note that the input arrays must be accessible when updating moments at
  each time step, so the array of inputs is declared static 
  in \ref cs_user_time_moments. Fo more complex inputs, we would have
  an array of inputs here; in this simple case, we could pass a simple
  call id as the input, casting from point to integer.

  \snippet cs_user_parameters-time_moments.c tmom_velocity_rotation

*/
