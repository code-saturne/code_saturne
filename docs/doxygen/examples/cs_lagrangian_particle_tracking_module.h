/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

  \page cs_lagrangian_particle_tracking_module Lagrangian volume statistical variable (uslaen.f90)
  
  \brief For the writing of the listing and the post-processing:
Average of the Lagrangian volume statistical variables.
Possible intervention of the user.

  
   \section cs_user_lagrangian_base_loc_var_lag Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_lagrangian_particle_tracking_module.f90 loc_var_dec

  \section cs_user_lagrangian_init_and_final Initialization and finalization

   The following initialization block needs to be added for the following examples:

  \snippet cs_user_lagrangian_particle_tracking_module.f90 allocate

  At the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_lagrangian_particle_tracking_module.f90 deallocate
 
  In theory Fortran 95 deallocates locally-allocated arrays automatically,
  but deallocating arrays in a symetric manner to their allocation is good pratice,
  and avoids using a different logic for C and Fortran.
  
  \section cs_user_lagrangian_standard_statistic Zone of standard statistics

 Pinpoint the cells where stats are to be calculated.

 \snippet cs_user_lagrangian_particle_tracking_module.f90 pinpoint_lag

 
General case:

 - Component X of the particle velocity: ivarl=ilvx \n
 - Component Y of the particle velocity: ivarl=ilvy \n 
 - Component Z of the particle velocity: ivarl=ilvz \n
 - Particle temperature: ivarl=iltp \n
 - Particle diameter: ivarl=ildp \n
 - Particle mass: ivarl= ilmp \n
 - Temperature of the coal particles: ivarl=ilhp(ilayer) \n
 - Mass of moisture of the coal particles: ivarl= ilmwat\n
 - Mass of reactive coal of the coal particles: ivarl= ilmch\n
 - Mass of coke of the coal particles: ivarl=ilmck\n
 - Diameter of the shrinking core of the coal particles: ivarl=ilmck\n
  except volume fraction (ivarl=ilfv) and sum of the statistical weights
  (ivarl=ilpd)


\snippet cs_user_lagrangian_particle_tracking_module.f90 general_lag

 Average
  \snippet cs_user_lagrangian_particle_tracking_module.f90 average_lag
Variance
  \snippet cs_user_lagrangian_particle_tracking_module.f90 variance_lag
<b>Volume fraction (ilfv)</b>
 \snippet cs_user_lagrangian_particle_tracking_module.f90 vol_lag
Average
  \snippet cs_user_lagrangian_particle_tracking_module.f90 average_lag_2
Variance
  \snippet cs_user_lagrangian_particle_tracking_module.f90 variance_lag_2
Sum of the statistical weights
  \snippet cs_user_lagrangian_particle_tracking_module.f90 sum_stat_weights
 

  \section cs_user_lagrangian_format_lag Format 
 
  \snippet cs_user_lagrangian_particle_tracking_module.f90 format_lag

  \section cs_user_lagrangian_user Zone of the user intervention: example 
Example 1: Statistic calculated in \ref uslast.90 and stored in the array statis.


  \snippet cs_user_lagrangian_particle_tracking_module.f90 example_1
  
*/
// __________________________________________________________________________________
