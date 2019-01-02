!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \file usthht.f90
!>
!> \brief Enthalpy-temperature conversion law definition.
!>
!>  Enthalpy    -> Temperature law (mode =  1) \n
!>  Temperature -> Enthalpy    law (mode = -1)
!>
!> See \subpage us_thht for examples.

subroutine usthht &
!================

 ( mode   , enthal , temper  )

!===============================================================================
!> \brief Enthalpy-temperature conversion law definition.
!>
!>  Enthalpy    -> Temperature law (mode =  1) \n
!>  Temperature -> Enthalpy    law (mode = -1)
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     mode          -1 : t -> h  ;   1 : h -> t
!> \param[in]     enthal        enthalpie
!> \param[in]     temper        temperature
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          mode

double precision enthal, temper

!===============================================================================


!     WARNING, the subroutine is called in loops:
!     =======                    ===============

!       Avoid parallel (collective) operations
!       =====

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

write(nfecra,9000)
call csexit (1)
!==========

!----
! Formats
!----

 9000 format(                                                           &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in enthalpy-temperature conversion',/,     &
'@    =======',/,                                                 &
'@',/,                                                            &
'@     The user subroutine ''usthht'' must be completed.',/,      &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!----
! End
!----

return
end subroutine usthht
