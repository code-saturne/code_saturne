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

/*-----------------------------------------------------------------------------*/

/*!

  \page cs_lagrangian_particle_tracking_module Parameters settings for lagrangian module

  \section cs_user_lagr_module_intro  Introduction

  This page gives some examples of settings for the stochastic lagrangian module.

  \section cs_user_lagr_module_h  Lagrangian module

  Particle tracking mode settings:

  \snippet cs_user_lagr_model.c particle_tracking_mode

  In case of restart

  \snippet cs_user_lagr_model.c particle_tracking_restart

  Specific models

  \snippet cs_user_lagr_model.c particle_tracking_specific_models

  Example of coal fouling

  \snippet cs_user_lagr_model.c coal_fouling_example

  Calculation features for the dispersed phases

  \snippet cs_user_lagr_model.c dispersed_phases

  Example of volume statistics

  \snippet cs_user_lagr_model.c V_statistics

  Options concerning the numerical treatment of the dispersed phase

  \snippet cs_user_lagr_model.c dispersed_phases_treatment

  Options concerning the treatment of specific forces

  \snippet cs_user_lagr_model.c specific_forces_treatment

  Brownian motion:

  \snippet cs_user_lagr_model.c Brownian_motion_activation

  Deposition model:

  \snippet cs_user_lagr_model.c deposition_model_activation

  Roughness resuspension model

  \snippet cs_user_lagr_model.c roughness_resuspension_model_activation

  Clogging model

  \snippet cs_user_lagr_model.c clogging_model_activation

  Deposit influence

  \snippet cs_user_lagr_model.c deposit_influence_activation

  Consolidation model:

  \snippet cs_user_lagr_model.c consolidation_model_activation

  Precipitation disolution model

  \snippet cs_user_lagr_model.c precipitation_disolution_model_activation

  Boundary statistics

  \snippet cs_user_lagr_model.c boundary_statistics



  \page cs_lagrangian_particle_tracking_bc User boundary condition definition for the Lagrangian model

  \section cs_user_lagr_boundary_conditions_h  Boundary conditions

  Lagrangian boundary conditions are based on boundary zone
  (\ref cs_boundary_zone_t) definitions. Additional information may be
  provided for Lagrangian boundary types and injections.

  As usual, definitions may be created using the GUI and extended
  with user functions.

  Access to the Lagrangian boundary conditions structure,
  which is necessary to most of the following examples, may be done as
  follows:

  \snippet cs_user_lagr_boundary_conditions.c lagr_bc_variables

  \subsection cs_user_lagr_boundary_conditions_h_zones  Boundary zones

  In this example, we assign rebound conditions to all boundary zones,
  except for an inlet and outlet type to specified zones.
  The type assigned is an integer based on the \ref cs_lagr_bc_type_t
  enumerator type.

  \snippet cs_user_lagr_boundary_conditions.c lagr_bc_define_type_1

  \subsection cs_user_lagr_boundary_conditions_h_injection  Injection sets

  In the following example, a first injection set for an inlet zone is defined.
  Note that newly injected particles may also be modified using the
  \ref cs_user_lagr_in function.

  \snippet cs_user_lagr_boundary_conditions.c lagr_bc_define_injection_1

  In the next example, a profile is assigned to the second injection set
  of an inlet zone (it is assumed this et was previously defined either
  through the GUI or user function).

  This requires first defining a profile definition function, matching
  the profile of \ref cs_lagr_injection_profile_compute_t.
  An example based on experimental profiles is given here:

  \snippet cs_user_lagr_boundary_conditions.c lagr_bc_profile_func_2

  Assigning the profile to the injection set simply requires
  assigning the function to the pointer in the injection set structure:

  \snippet cs_user_lagr_boundary_conditions.c lagr_bc_define_injection_2

  An optional user-defined input function may also be associated.

  \subsection cs_user_lagr_boundary_conditions_h_interactions  Boundary-particle interactions

  It is also possible to decide of the behavior of particle when they
  encounter a boundary (this boundary has to be of type \ref CS_LAGR_BC_USER).

  In the following example, the particle is simply deposited and marked for
  elimination:

  \snippet cs_user_lagr_boundary_conditions.c update



  \section cs_user_lagr_volume_conditions_h  Volume conditions

  Lagrangian volume conditions are based on volume zone
  (\ref cs_volume_zone_t) definitions. Additional information may be
  provided for Lagrangian injections.

  As usual, definitions may be created using the GUI and extended
  with user functions.

  Access to the Lagrangian volume conditions structure,
  which is necessary to most of the following examples, may be done as
  follows:

  \snippet cs_user_lagr_volume_conditions.c lagr_vol_cond_variables

  \subsection cs_user_lagr_volume_conditions_h_injection  Injection sets

  In the following example, we inject 1 particle set at each time step:
  \snippet cs_user_lagr_volume_conditions.c lagr_vol_define_injection_1

  \snippet cs_user_lagr_volume_conditions.c lagr_vol_define_injection_2
  In the following example, we inject 2 particle sets at computation
  initialization (i.e. at the first time step of a computation sequence
  in which the Lagrangian module is activated).
  Note that newly injected particles may also be modified using the
  \ref cs_user_lagr_in function.

  \snippet cs_user_lagr_volume_conditions.c lagr_vol_define_injection_2

  \subsection cs_user_lagr_volume_conditions_h_force  External force

  User definition of an external force field acting on the particles.

  It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
  By default gravity and drag force are the only forces acting on the particles
  (the gravity components gx gy gz are assigned in the GUI or in usipsu)

  \snippet cs_user_lagr_particle.c lagr_ef

  \subsection cs_user_lagr_volume_conditions_h_imposed_motion Impose motion on a particle

  Impose the motion of a particle flagged CS_LAGR_PART_IMPOSED_MOTION by modifying the particle position and its velocity.

  \snippet cs_user_lagr_particle.c lagr_imposed_motion

  \subsection cs_user_lagr_volume_conditions_h_inj_particles Modification of newly injected particles

  User modification of newly injected particles.

  This function is called after the initialization of the new particles in
  order to modify them according to new particle profiles (injection
  profiles, position of the injection point, statistical weights,
  correction of the diameter if the standard-deviation option is activated).

  This function is called for each injection zone and set. Particles
  with ids between  pset->n_particles and  n_elts are initialized
  but may be modidied by this function.

  \snippet cs_user_lagr_particle.c lagr_inj

  Here is another example of the modification of newly injected particles:

  \snippet cs_user_lagr_particle-coal.c lagr_inj_example_coal

  \page cs_lagrangian_particle_tracking_physical_properties Particle relaxation time and thermal relaxtion for the Lagrangian model

  \section cs_user_lagr_module_time_relaxation Calculation of the particle relaxation time

  Modification of the calculation of the particle relaxation time
  with respect to the chosen formulation for the drag coefficient

  This function is called in a loop on the particles, so be careful to avoid too costly operations.

  \tau_c = \frac{m_p{C_p}_p}{PId_p^2h_e}

  \tau_c  : Thermal relaxation time (value to be computed)

  m_p    : Particle mass

  {C_p}_p   : Particle specific heat

  d_p    : Particle diameter

  h_e    : Coefficient of thermal exchange

  The coefficient of thermal exchange is calculated from a Nusselt number,
  itself evaluated by a correlation (Ranz-Marshall by default)

  \nu =  \frac{h_ed_p}{\lambda} = 2 + 0.55{\Re_e}_p^{0.5}P_{rt}^{0.33}

  \lambda : Thermal conductivity of the carrier field

  {\Re_e}_p     : Particle Reynolds number

  P_{rt}    : Prandtl number


  In the next example we compute the relaxation time with two different formulation of the drag coefficient:

  \snippet cs_user_lagr_particle.c lagr_particle_relax_time

  \section cs_user_lagr_module_thermal_relaxation Computation of the thermal relaxation time of the particles

  Modification of the computation of the thermal relaxation time
  of the particles with respect to the chosen formulation of
  the Nusselt number.

  This function is called in a loop on the particles, so be careful to avoid too costly operations.

  \snippet cs_user_lagr_particle.c lagr_thermal_relax_time



  \page cs_user_lagr_extra_operations User extra operations for the Lagrangian module

  \section cs_user_lagr_extra_operations_intro Introduction

  This page provides an example that may be used or adapted to perform extra or advanced extra-operations within the Lagrangian module.

  \section cs_user_lagr_extra_operations_example Example

  First we initialize some variables:

  \snippet cs_user_lagr_particle.c lagr_init

  In the next example we compute the particle mass flow rate on 4 planes

  \snippet cs_user_lagr_particle.c lagr_example



  \page cs_user_lagr_sde SDE within the Lagrangian model

  \section cs_user_lagr_sde_intro Introduction

  User integration of the SDE for the user-defined variables.

  The variables are constant by default. The SDE must be of the form:

  \f[
  \frac{dT}{dt}=\frac{T - PIP}{Tca}
  \f]

  T:   particle attribute representing the variable

  Tca: characteristic time for the sde to be prescribed in the array auxl1

  PIP: coefficient of the SDE (pseudo RHS) to be prescribed in the array auxl2.

  If the chosen scheme is first order (nordre=1) then, at the first
  and only call pip is expressed as a function of the quantities of
  the previous time step (contained in the particle data).

  If the chosen scheme is second order (nordre=2)
  then, at the first call (nor=1) pip is expressed as a function of
  the quantities of the previous time step, and at the second passage
  (nor=2) pip is expressed as a function of the quantities of the
  current time step.

  \section cs_user_lagr_sde_example Example

  Example of the integration of the SDE (Stochastic Differential Equation).

  \snippet cs_user_lagr_particle.c lagr_SDE


*/
