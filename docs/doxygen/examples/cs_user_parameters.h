/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*!
  \page parameters Input of calculation parameters (C functions in cs_user_parameters.c)

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

  \section cs_user_parameters_h_cs_user_parameters General options (\ref cs_user_parameters)

  Most definitions should be done in the \ref cs_user_parameters
  function (Fortran equivalent is \ref usipsu subroutine).

  Choose a turbulent model among the available models

  \snippet cs_user_parameters-base.c turbulence_model_choice

  Coupled solver for Rij components (when iturb=30, 31 or 32): 0 to switch off, 1 to switch on

  \snippet cs_user_parameters-base.c Rij_coupled_solver_choice

  To set a thermal model (0: none, 1: temperature, 2: entyhalpy, 3: total energy)

  \snippet cs_user_parameters-base.c thermal_model_choice

  Volume of Fluid model with mass transfer Merkle model to take into account
  vaporization / condensation

  \snippet cs_user_parameters-base.c enable_cavit

  To activate ALE (Arbitrary Lagrangian Eulerian) method
  (CS_ALE_NONE: switch off, CS_ALE_LEGACY: legacy solver, CS_ALE_CDO: CDO solver)

  \snippet cs_user_parameters-base.c ALE_activation

  The user can add a scalar to be solved

  \snippet cs_user_parameters-base.c scalars_addition

  After adding a scalar, the user can add the variance of this scalar

  \snippet cs_user_parameters-base.c scalars_variance_addition

  Add a user property defined on a mesh location (cells, interior faces, boundary faces or vertices).

  \snippet cs_user_parameters-base.c user_property_addition

  Choose a time step option

  \snippet cs_user_parameters-base.c time_stepping_options

  Choose a reference time step

  \snippet cs_user_parameters-base.c ref_time_step

  To set a duration

  \snippet cs_user_parameters-base.c duration

  For example, to change the log (run_solver.log) verbosity of all the variables:

  \snippet cs_user_parameters-base.c param_log_verbosity

  For example, to change limiters for the convective scheme for a given scalar:

  \snippet cs_user_parameters-base.c param_var_limiter_choice

  One can also choose the percentage of upwind blending when using the slope test

  \snippet cs_user_parameters-base.c param_var_blend_st

  If one wants to declare a scalar as buoyant (i.e. the density depends on this scalar
  through a given equation of state) and add it in the velocity-pressure optional
  sub-iterations (PISO), one can set the dedicated keyword:

  \snippet cs_user_parameters-base.c param_var_is_buoyant

  If one wants to activate drift terms on a transported scalar:

  \snippet cs_user_parameters-base.c param_var_drift


  To set options for the solving of the Stokes problem (velocity pressure coupling
  see \ref cs_stokes_model_t structure):

  \snippet cs_user_parameters-base.c param_stokes_model

  To set options for the PISO-like sub iterations (velocity pressure coupling
  see \ref cs_piso_t structure):

  \snippet cs_user_parameters-base.c param_piso


  \subsection cs_user_parameters_h_cs_user_parameters1 Special fields and activate some automatic post-processings

  For example, to force the presence of a boundary temperature field, which
  may be useful for postprocessing:

  \snippet cs_user_parameters-base.c param_force_b_temperature

  To add boundary values for all scalars, the
  following code can be added:

  \snippet cs_user_parameters-base.c param_var_boundary_vals_1

  Enforce existence of 'tplus' and 'tstar' fields, so that
  a Nusselt number may be computed using the
  \ref post_boundary_nusselt subroutine.
  When postprocessing this quantity is activated, those fields
  are present, but if we need to compute them in the
  \ref cs_user_extra_operations user subroutine without postprocessing them,
  forcing the definition of these fields to save the values computed
  for the boundary layer is necessary.

  \snippet cs_user_parameters-base.c param_force_yplus

  You can activate the post-processing of the Q-criterion on the whole domain
  mesh with:

  \snippet cs_user_parameters-base.c param_var_q_criterion

  Save contribution of slope test for variables in special fields.
  These fields are automatically created, with postprocessing output enabled,
  if the matching variable is convected, does not use a pure upwind scheme,
  and has a slope test (the slope_test_upwind_id key value for a given
  variable's field is automatically set to the matching postprocessing field's
  id, or -1 if not applicable).

  \snippet cs_user_parameters-base.c param_post_slop_test

  You can activate the post-processing of clipping on turbulent quantities
  on the whole domain mesh with:

  \snippet cs_user_parameters-base.c param_var_rij_clipping


  \section cs_user_parameters_h_finalize_setup Input-output related examples (usipes)

  \subsection cs_user_parameters_h_example_base Basic options

  Frequency of log output.

  \snippet cs_user_parameters-base.c setup_log

  Change a property's label (here for density, first checking if it
  is variable). A field's name cannot be changed, but its label,
  used for logging and postprocessing output, may be redefined.

  \snippet cs_user_parameters-base.c setup_label

  \subsection cs_user_parameters_h_example_post Postprocessing output

  Activate or deactivate probes output.

  \snippet cs_user_parameters-base.c setup_post

  Probes for Radiative Transfer (Luminance and radiative density flux vector).

  \snippet cs_user_parameters-base.c setup_post_lum


  \section cs_user_parameters_h_cs_user_model_cdo Base model for CDO/HHO schemes

  \subsection cs_user_parameters_h_cdo_activation Activation of CDO/HHO schemes

  Two modes are available \ref CS_DOMAIN_CDO_MODE_ONLY or \ref CS_DOMAIN_CDO_MODE_WITH_FV

  CDO/HHO schemes can be activated within this function as follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_activation

  \subsection cs_user_parameters_h_cdo_domain_boundary Domain boundary with CDO/HHO schemes

  Several types of domain boundaries can be defined. There are gathered in \ref cs_domain_boundary_type_t
  The definition of the domain boundaries for CDO/HHO schemes can be specified
  as follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_domain_boundary

  \subsection cs_user_parameters_h_cdo_domain_output Output with CDO/HHO schemes

  The management of the level and frequency of details written by the code
  can be specified for CDO/HHO schemes as follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_domain_output

  \subsection cs_user_parameters_h_cdo_time_step Time step with CDO/HHO schemes

  The management of the time step with CDO/HHO schemes can be specified as
  follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_time_step

  \subsection cs_user_parameters_h_cdo_walldistance Wall distance with CDO/HHO schemes

  The computation of the wall distance with CDO schemes is performed as follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_wall_distance

  \subsection cs_user_parameters_h_cdo_add_user_eq Add a user-defined equation with CDO/HHO schemes

  The add of a user-defined equation solved by CDO/HHO schemes is specified as
  follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_equation

  \subsection cs_user_parameters_h_cdo_add_user_pty Add user-defined properties with CDO/HHO schemes

  The add of a new user-defined property with CDO/HHO schemes is specified as follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_properties

  If you want to compute the Fourier number related to a given property in an unsteady simulation
  \snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_properties_opt

  \subsection cs_user_parameters_h_cdo_add_user_adv_field Add user-defined advection field with CDO/HHO schemes

  The definition of an advection field allows one to handle flows with a frozen velocity field or
  the transport of scalar quantities without solving the Navier-Stokes system.
  The add of a new user-defined advection field with CDO/HHO schemes is specified as follows:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_adv_field

  If you need to activate options related to advection fields, you can also specify

  \snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_adv_field_opt

  \subsection cs_user_parameters_h_cdo_gwf Settings related to the groundwater flow module with CDO/HHO schemes

  Activation of the groundwater flow module. The second argument is either 0 if no option is needed or a list of flags among
  \ref CS_GWF_GRAVITATION, \ref CS_GWF_RICHARDS_UNSTEADY, \ref CS_GWF_SOIL_PROPERTY_UNSTEADY, \ref CS_GWF_SOIL_ALL_SATURATED

  \snippet cs_user_parameters-cdo-gwf.c param_cdo_activate_gwf

  Add soils (settings of the soil is performed in \ref cs_user_gwf_setup). Soils have to be added before adding tracers.

  \snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_add_soil

  Add tracer equations which correspond to a transport equation using
  the darcean flux as the advection field. This call implies the
  creation of a new equation related to a new variable field name
  given as arguments. Advanced tracer equations can be defined using
  \ref cs_gwf_add_tracer_user.

  \snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_add_tracer

  \section cs_user_parameters_h_cs_user_parameters Define or modify general numerical and physical user parameters

  \subsection cs_user_parameters_h_cdo_param_eq CDO/HHO schemes

  Modifiy the numerical parameters related to a given equation.

  \snippet cs_user_parameters-cdo-condif.c param_cdo_numerics

  \section cs_user_parameters_h_cs_user_cdo_finalize_setup Finalize the set-up for CDO/HHO schemes

  \subsection cs_user_parameters_h_cdo_set_pty Set up properties with CDO/HHO schemes

  When a property has been added, the second step is to define this property. According to the type
  of property (isotropic, orthotropic or anisotropic) definitions differ.
  Here are two examples:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_setup_property

  \subsection cs_user_parameters_h_cdo_set_advfield Set up user-defined advection field with CDO/HHO schemes

  When an advection field has been added, the second step is to define this advection field.
  Here are is an example of definition using an anlytic function and the activation of optional features:

  \snippet cs_user_parameters-cdo-condif.c param_cdo_setup_advfield

  \subsection cs_user_parameters_h_cdo_set_bcs Set up the boundary conditions for an equation

  \snippet cs_user_parameters-cdo-condif.c param_cdo_setup_bcs

  \subsection cs_user_parameters_h_cdo_add_eq_terms Add terms to an equation

  Add terms like diffusion term, advection term, unsteady term, reaction terms or source terms.

  \snippet cs_user_parameters-cdo-condif.c param_cdo_add_terms

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

  \subsection cs_user_parameters_h_sles_verbosity_1 Changing the verbosity

  By default, a linear solver uses the same verbosity as its matching variable,
  and is not verbose for non-variable quantities. The verbosity
  may be specifically set for linear system resolution, as shown in
  the following example:

  \snippet cs_user_parameters-linear_solvers.c sles_verbosity_1

  \subsection cs_user_parameters_h_sles_viz_1 Example: error visualization

  The following example shows how to activate local error visualization
  output (here for velocity and pressure).

  \snippet cs_user_parameters-linear_solvers.c sles_viz_1

  \subsection cs_user_parameters_h_sles_mgp_1 Example: advanced multigrid settings

  The following example shows how to set advanced settings for the
  multigrid solver used for the pressure solution.

  \snippet cs_user_parameters-linear_solvers.c sles_mgp_1

  \subsection cs_user_parameters_h_sles_mgp_2 Example: multigrid preconditioner settings

  The following example shows how to use multigrid as a preconditioner for a
  conjugate gradient solver (for the pressure correction), and set advanced
  settings for that multigrid preconditioner.

  \snippet cs_user_parameters-linear_solvers.c sles_mgp_2

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

  \subsection cs_user_parameters_h_sles_plot Plotting solver convergence

  The following example shows how to activate convergence plotting
  for built-in iterative or multigrid solvers.

  \snippet cs_user_parameters-linear_solvers.c sles_plot_1

  Plots will appear as CSV (comma-separated value) files in the
  monitoring subdirectory.

  \subsection cs_user_parameters_h_sles_petsc Using PETSc

  The following example shows how to setup a solver to use the PETSc
  library, if the code was built with PETSc support.

  General options (those passed to PETSc through command line options)
  may be defined directly in \ref cs_user_linear_solvers, for example:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_1

  A specific system may be set up to use PETsc, as is shown
  here for the pressure variable:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_2

  The basic matrix format to be used by PETSc is defined at this stage, using
  a PETSc MatType string (see PETSc documentation).
  Further options may be defined in a setup hook function, as follows:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_hook_1

  If no additional settings are required, the matching parameter
  in \ref cs_sles_petsc_define may be set to NULL.

  \subsubsection cs_user_parameters_h_sles_petsc_gamg Using PETSc with GAMG

  To use PETSc's GAMG (geometric-algebraic multigrid) preconditioner,
  the following general options should be set:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_gamg_1

  Setting GAMG-preconditioned PCG for the pressure may be done as in
  the previous option, ensuring here that a matrix structure
  native to PETSc is used (SEQAIJ, MPIAIJ, or some other AIJ variant),
  so all required matrix operations are available:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_gamg_2

  With the associated setup hook function:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_hook_gamg

  \subsubsection cs_user_parameters_h_sles_petsc_bamg Using PETSc with HYPRE

  To use HYPRE's Boomer AMG as a PETSc preconditioner,
  the following general options should be set:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_bamg_1

  Setting BoomerAMG-preconditioned PCG for the pressure may be done as in
  the previous option, ensuring here that a matrix structure
  native to PETSc is used (SEQAIJ, MPIAIJ, or some other AIJ variant),
  so all required matrix operations are available:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_bamg_2

  With the associated setup hook function:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_hook_bamg

  \subsubsection cs_user_parameters_h_sles_petsc_add Additional PETSc features

  Many additional features are possible with PETSc; for example,
  the following setup hook also outputs a view of the matrix, depending
  on an environment variable, \c CS_USER_PETSC_MAT_VIEW, which may
  take values \c DEFAULT, \c DRAW_WORLD, or \c DRAW:

  \snippet cs_user_parameters-linear_solvers.c sles_petsc_hook_view

  \section cs_user_parameters_h_cs_user_moments  Time moment related options

  Code_Saturne allows the calculation of temporal means or variances,
  either of expressions evaluated through a user function, or
  of expressions of the type \f$<f_1*f_2...*f_n>\f$. The variables
  may be fields or field components. This is done calling either
  through the GUI, or in the user function \ref cs_user_time_moments.
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

  To count defined moments, use the \ref cs_time_moment_n_moments function,
  whether from Fortran or C. To access the matching fields, use
  \ref cs_c_bindings::time_moment_field_id "time_moment_field_id" in Fortran,
  or \ref cs_time_moment_get_field in C.

  \section cs_user_parameters_h_examples Examples

  \subsection cs_user_parameters_h_example_1 Example 1

  In the following example, we define a moment for the mean velocity.
  All components are used (component -1 means all components),
  so the moment is a vector.

  \snippet cs_user_parameters-time_moments.c tmom_u

  In the next example, we define the variance of the vector velocity.
  All components are used again (component -1 means all components),
  so the moment is a tensor.

  \snippet cs_user_parameters-time_moments.c tmom_variance_u

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

  In \ref cs_user_time_moments, we can now assign that function to a moments
  definition:

  \snippet cs_user_parameters-time_moments.c tmom_simple_sum

  \subsection cs_user_parameters_h_example_5 Example 5

  This next example illustrates the use of another user-defined function
  to evaluate expressions. Here, we compute the moment of the thermal flux
  at the boundary. We also that we compute both means and variances here.

  \snippet cs_user_parameters-time_moments.c tmom_b_thermal_flux_data

  In \ref cs_user_time_moments, we assign that function to a moments
  definition:

  \snippet cs_user_parameters-time_moments.c tmom_b_thermal_flux

  \subsection cs_user_parameters_h_example_6 Example 6

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

  To activate means for all variables:
  \snippet cs_user_parameters-time_moments.c tmom_all_variables

  \section cs_user_parameters_h_cs_user_fans  Fan modelling options

  Code_Saturne allows modelling of some circular fans as volume
  regions, defined by simple geometric characteristics, and modeled
  as explicit momentum source terms in those regions.

  Fan pressure characteristic curves are defined as a 2nd order
  polynomial, and a torque may also be specified. For correct results,
  it is important that the mesh match the fan dimensions and placement
  (thickness, hub, blades, and total radius).

  The following example shows how a fan may be defined:

  \snippet cs_user_parameters-fans.c fan_user_1

*/
