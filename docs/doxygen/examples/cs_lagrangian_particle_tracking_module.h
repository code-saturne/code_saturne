/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

  \page cs_lagrangian_particle_tracking_module Lagrangian module

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

  In the following example, we inject 2 particle sets at computation
  initialization (i.e. at the first time step of a computation sequence
  in which the Lagrangian module is activated).
  Note that newly injected particles may also be modified using the
  \ref cs_user_lagr_in function.

  \snippet cs_user_lagr_volume_conditions.c lagr_vol_define_injection_1
*/
