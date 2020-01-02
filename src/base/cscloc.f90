!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file cscloc.f90
!> \brief Coupling interfaces localization (with FVM).
!>
!------------------------------------------------------------------------------

subroutine cscloc

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use cplsat

!===============================================================================

implicit none

! Arguments

! Local variables

integer          numcpl

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

ipass  = ipass + 1

do numcpl = 1, nbrcpl

  ! On localise au premier passage ou si l'un des maillages
  ! du couplage est mobile ou avec ALE (cf CSCINI).
  if (ipass.eq.1.or.imajcp(numcpl).eq.1) then

    ! Localisation proprement dite
    call  defloc ( numcpl )
    !===========

  endif

enddo


!--------
! FORMAT
!--------

!----
! FIN
!----

return
end subroutine
