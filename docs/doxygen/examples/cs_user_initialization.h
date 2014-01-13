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
  \page cs_user_initialization cs_user_initialization.f90


  \section intro Introduction

  This page provides several examples of code blocks that may be used
  to initialize variables in \ref cs_user_initialization.

  Theese subroutines are called at beginning of the computation
  (restart or not) before the loop time step.

  Theese subroutines enable to initialize or modify (for restart)
  unkown variables and time step values.

  rom and viscl values are equal to ro0 and viscl0 or initialize
  by reading the restart file:
  - viscls and cp variables (when there are defined) have no value
  - excepted if they are read from a restart file.

  Physical quantities are defined in the following arrays:
   - propce (physical quantities defined at cell center),
   - propfb (physical quantities defined at border face center).
 
  Examples:
   - propce(iel, ipproc(irom  )) means rom  (iel)
   - propce(iel, ipproc(iviscl)) means viscl(iel)
   - propce(iel, ipproc(icp   )) means cp   (iel)
   - propce(iel, ipproc(ivisls(iscal))) means visls(iel, iscal)
   - propfb(ifac, ipprob(irom )) means romb  (ifac)

  Modification of the behaviour law of physical quantities (rom, viscl,
  viscls, cp) is not done here. It is the purpose of the user subroutine
  usphyv.

 Cells identification

  Cells may be identified using the 'getcel' subroutine.
  The syntax of this subroutine is described in the
  \ref cs_user_boundary_conditions subroutine,
  but a more thorough description can be found in the user guide.

  \section cs_user_bc_examples Initialization examples
  Here is the list of examples dedicated to different physics:

   - \subpage base
   - \subpage atmospheric
   - \subpage compressible
   - \subpage cooling_towers
   - \subpage electric_arcs
   - \subpage fuel
   - \subpage gas_3ptchem
   - \subpage gas_ebu
   - \subpage gas_libby_williams
   - \subpage pulverized_coal
   - \subpage time_step
   - \subpage unified_combustion

*/
// __________________________________________________________________________________
/*!

  \page base Basic example

  \section base Basic example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-base.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays must be allocated.

  \snippet cs_user_initialization-base.f90 alloc

  \subsection s_init Initialization

   isca(1) is the number related to the first user-defined scalar variable.
   rtp(iel,isca(1)) is the value of this variable in cell number iel.

   ONLY done if there is no restart computation.

  \snippet cs_user_initialization-base.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_initialization-base.f90 finalize

*/
// __________________________________________________________________________________
/*!

  \page atmospheric Atmospheric example

  \section atmospheric Atmospheric example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-atmospheric.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-atmospheric.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page compressible Compressible example

  \section compressible Compressible example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-compressible.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays must be allocated.

  \snippet cs_user_initialization-base.f90 alloc

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-compressible.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work arrays:

  \snippet cs_user_initialization-compressible.f90 finalize

*/
// __________________________________________________________________________________
/*!

  \page cooling_towers cooling_towers example

  \section cooling_towers cooling_towers example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-cooling_towers.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-cooling_towers.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page electric_arcs Electric arcs example

  \section electric_arcs Electric arcs example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-electric_arcs.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  Some controls are made during initialization:   

  \snippet cs_user_initialization-electric_arcs.f90 init1

  Classical initialization:

  \snippet cs_user_initialization-electric_arcs.f90 init2

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page fuel Fuel example

  \section fuel Fuel example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-fuel.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-fuel.f90 init

  \subsection finalize Finalization

  There is no work array in this subroutine, thus nothing to do.

*/
// __________________________________________________________________________________
/*!

  \page gas_3ptchem Gas 3 PTCHEM example

  \section gas_3ptchem Gas 3 PTCHEM example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-gas_3ptchem.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-gas_3ptchem.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page gas_ebu Gas EBU example

  \section gas_ebu Gas EBU example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-gas_ebu.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-gas_ebu.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page gas_libby_williams Libby-Williams gas example

  \section gas_libby_williams  Libby-Williams gas example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-gas_libby_williams.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-gas_libby_williams.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page pulverized_coal pulverized_coal example

  \section pulverized_coal pulverized_coal example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-pulverized_coal.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-pulverized_coal.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page time_step Time step modification

  \section time_step Time step modification

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-time_step.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-time_step.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
// __________________________________________________________________________________
/*!

  \page unified_combustion Unified combustion coal example
  isca(1) is the number related to the first user-defined scalar variable.
  rtp(iel,isca(1)) is the value of this variable in cell number iel.
  ONLY done if there is no restart computation

  \section unified_combustion Unified combustion coal example

  \subsection loc_var_dec Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_initialization-unified_combustion_coal.f90 loc_var_dec

  \subsection alloc Allocation

  Before user initialization, work arrays lstelt must be allocated, 
  like in basic example.

  \subsection s_init Initialization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_initialization-unified_combustion_coal.f90 init

  \subsection finalize Finalization

  At the end of the subroutine, it is recommended to deallocate the work array lstelt,
  like in basic example.

*/
