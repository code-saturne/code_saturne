!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_parameters-cavitation.f90
!>
!> \brief Cavitation parameters example.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------

!===============================================================================

subroutine usipsu &
!================

 ( nmodpp )


!===============================================================================
! Purpose:
! -------

! User subroutine for the input of additional user parameters.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nmodpp           ! i  ! <-- ! number of active specific physics models       !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, jj, imom, kscmin, kscmax

!===============================================================================

!     This subroutine allows setting parameters

!       which do not already appear in the other subroutines of this file.


!     It is possible to add or remove parameters.


!     The number of physical properties and variables is known here.


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

!-------------------------------------------------------------------------------
!>
!> \file cs_user_parameters-cavitation.f90
!>
!>
!>   Cavitation module:
!>   ------------------
!>
!>      Model based on a homogeneous mixture. The physical properties (density
!>        and dynamic viscosity) of the mixture depends on a resolved void
!>        fraction and constant reference properties of the liquid phase and
!>        the gas phase.
!>
!>    User example subroutine for cavitation module parameter definitions
!>
!-------------------------------------------------------------------------------


!===============================================================================
! 1. Definition of the homogeneous mixture physical properties
!===============================================================================

! --- Reference density, in kg/m3, and molecular viscosity, kg/(m s), of the
!       liquid phase

!< [phprop_l]
rol = 1.d3
mul = 1.d-3
!< [phprop_l]

! --- Reference density, in kg/m3, and molecular viscosity, kg/(m s), of the
!       gas phase

!< [phprop_g]
rov = 1.d0
muv = 1.d-5
!< [phprop_g]

!===============================================================================
! 2. Model parameters of the vaporization term (Merkle model)
!===============================================================================

! ---  Reference saturation pressure in kg/(m s2)

!< [presat]
presat = 2.d3
!< [presat]

! ---  Reference length, in meters, and velocity scales, in m/s, of the flow

!< [scales_inf]
linf = 0.1d0
uinf = 1.d0
!< [scales_inf]

!===============================================================================
! 3. Interaction with turbulence
!===============================================================================

! --- Eddy-viscosity correction (Reboud et al correction)
!      0: deactivated
!      1: activated

!< [reboud_activ]
icvevm = 1
!< [reboud_activ]

!----
! End
!----

end subroutine usipsu
