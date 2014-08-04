/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

  \section intro Introduction

  C user functions for definition of model options and calculation parameters.
    These subroutines are called in all cases.

  If the Code_Saturne GUI is used, this file is not required (but may be
    used to override parameters entered through the GUI, and to set
    parameters not accessible through the GUI).

  Several functions are present in the file, each destined to defined
    specific parameters.


  \section cs_user_model  Base model related options

  Definition of user variables or properties should be defined here,
  if not already done throught the GUI.

  \section cs_user_moments  Time moment related options

  Code_Saturne allows the calculation of temporal means or variances,
  either of expressions evaluated through a user function, or
  of expressions of the type \f$f_1*f_2...*f_n\f$. The variables
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

  \section examples Examples

  \subsection example_1 Example 1

  In the following example, we define a moment for the mean velocity.
  All components are used (component -1 means all components),
  so the moment is a vector.

  \snippet cs_user_parameters-time_moments.c tmom_u

  \subsection example_2 Example 2

  In the next example, we multiply the expression by the density.
  As the density is of dimension 1, and the velocity of dimension 3,
  the resulting moment is of dimension 3.

  \snippet cs_user_parameters-time_moments.c tmom_rho_u

  \subsection example_3 Example 3

  In the next example, we define a product of several field components,
  all of dimension 1, as we consider only the x and y components of the
  velocity; for the density, we cas use either component 0 or -1 (all),
  since the field is scalar.

  This moment's computation is also restarted at each time step.

  \snippet cs_user_parameters-time_moments.c tmom_rho_u_v

  \subsection example_4 Example 4

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

  \subsection example_5 Example 5

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
