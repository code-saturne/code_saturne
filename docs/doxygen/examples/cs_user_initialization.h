/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
  \page user_initialization cs_user_initialization.f90


  \section user_initialization_intro Introduction

  This page provides several examples of code blocks that may be used
  to initialize variables in \ref cs_user_initialization.

  These subroutines are called at beginning of the computation
  (restart or not) before the loop time step.

  These subroutines enable to initialize or modify (for restart)
  unkown variables and time step values.

  rom and viscl values are equal to ro0 and viscl0 or initialize
  by reading the restart file:
  - viscls and cp variables (when there are defined) have no value
  - excepted if they are read from a restart file.

  Modification of the behaviour law of physical quantities (rom, viscl,
  viscls, cp) is not done here. It is the purpose of the user subroutine
  usphyv.

 Cells identification

  Cells may be identified using the 'getcel' subroutine.
  The syntax of this subroutine is described in the
  \ref cs_user_boundary_conditions subroutine,
  but a more thorough description can be found in the user guide.

  \section cs_user_init_examples Initialization examples
  Here is the list of examples dedicated to different physics:

   - \subpage user_initialization_base
   - \subpage user_initialization_atmospheric
   - \subpage user_initialization_compressible
   - \subpage user_initialization_electric_arcs
   - \subpage user_initialization_fuel
   - \subpage user_initialization_gas_3ptchem
   - \subpage user_initialization_gas_ebu
   - \subpage user_initialization_gas_libby_williams
   - \subpage user_initialization_pulverized_coal
   - \subpage user_initialization_time_step
   - \subpage user_initialization_unified_combustion

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_base Basic example

  \section base Basic example

  \subsection user_initialization_base_loc_var_dec Local variables to be added

  \subsection user_initialization_base_s_init Initialization

  One can get any field using \ref cs_field_by_name function (use
  \ref cs_field_by_name_try if one is not sure the field exists).
  "scalar1" is the name related to the first user-defined scalar variable.
  \c f->val[\c cell_id] is the value of this variable in cell number \c cell_id.

  ONLY done if there is no restart computation.

  \snippet cs_user_initialization-base.c init

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_atmospheric Atmospheric example

  \section atmospheric Atmospheric example

  \subsection user_initialization_atmo_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-atmospheric.f90 loc_var_dec

  \subsection user_initialization_atmo_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_atmo_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-atmospheric.f90 init

  \subsection user_initialization_atmo_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_compressible Compressible example

  \section compressible Compressible example

  \subsection user_initialization_comp_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-compressible.f90 loc_var_dec

  \subsection user_initialization_comp_alloc Allocation

  Before user initialization, work arrays must be allocated.

  \snippet cs_user_initialization-base.f90 alloc

  \subsection user_initialization_comp_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-compressible.f90 init

  \subsection user_initialization_comp_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work arrays:

  \snippet cs_user_initialization-compressible.f90 finalize

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_electric_arcs Electric arcs example

  \section electric_arcs Electric arcs example

  \subsection user_initialization_ea_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-electric_arcs.f90 loc_var_dec

  \subsection user_initialization_ea_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_ea_s_init Initialization

  Classical initialization:

  \snippet cs_user_initialization-electric_arcs.f90 init2

  \subsection user_initialization_ea_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_fuel Fuel example

  \section fuel Fuel example

  \subsection user_initialization_fuel_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-fuel.f90 loc_var_dec

  \subsection user_initialization_fuel_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_fuel_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-fuel.f90 init

  \subsection user_initialization_fuel_finalize Finalization

  There is no work array in this subroutine, thus nothing to do.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_gas_3ptchem Gas 3 PTCHEM example

  \section gas_3ptchem Gas 3 PTCHEM example

  \subsection user_initialization_gas_3p_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-gas_3ptchem.f90 loc_var_dec

  \subsection user_initialization_gas_3p_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_gas_3p_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-gas_3ptchem.f90 init

  \subsection user_initialization_gas_3p_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_gas_ebu Gas EBU example

  \section gas_ebu Gas EBU example

  \subsection user_initialization_gas_ebu_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-gas_ebu.f90 loc_var_dec

  \subsection user_initialization_gas_ebu_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_gas_ebu_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-gas_ebu.f90 init

  \subsection user_initialization_gas_ebu_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_gas_libby_williams Libby-Williams gas example

  \section gas_libby_williams  Libby-Williams gas example

  \subsection user_initialization_gas_lb_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-gas_libby_williams.f90 loc_var_dec

  \subsection user_initialization_gas_lb_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_gas_lb_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-gas_libby_williams.f90 init

  \subsection user_initialization_gas_lb_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_pulverized_coal pulverized_coal example

  \section pulverized_coal pulverized_coal example

  \subsection user_initialization_coal_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-pulverized_coal.f90 loc_var_dec

  \subsection user_initialization_coal_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_coal_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-pulverized_coal.f90 init

  \subsection user_initialization_coal_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_time_step Time step modification

  \section time_step Time step modification

  \subsection user_initialization_time_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-time_step.f90 loc_var_dec

  \subsection user_initialization_time_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_time_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-time_step.f90 init

  \subsection user_initialization_time_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page user_initialization_unified_combustion Unified combustion coal example
  One can get any field using \ref field_get_val_s_by_name function.
  \c cvar_*(iel) is the value of this variable in cell number \c iel.
  ONLY done if there is no restart computation

  \section unified_combustion Unified combustion coal example

  \subsection user_initialization_comb_loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-unified_combustion_coal.f90 loc_var_dec

  \subsection user_initialization_comb_alloc Allocation

  Before user initialization, work arrays lstelt must be allocated,
  like in basic example.

  \subsection user_initialization_comb_s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-unified_combustion_coal.f90 init

  \subsection user_initialization_comb_finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
