!-------------------------------------------------------------------------------

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

subroutine csopli &
!================

 (infecr, isuppr, ierror)

!===============================================================================
! Purpose:
! -------

!    Open log files using Fortran IO.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! infecr           ! i  ! <-- ! value to assign to nfecra                      !
! isuppr           ! i  ! <-- ! supress output if ~                            !
! ierror           ! i  ! --> ! error code                                     !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor

!===============================================================================

implicit none

! Arguments

integer infecr, isuppr, ierror

! Local variables

character(len=64) :: name

!===============================================================================

ierror = 0

nfecra = infecr

if (nfecra .eq. 6) return

call cslogname(len(name), name)

if (isuppr .eq. 0) then
  open(file=name, unit=nfecra, form='formatted', status='old',   &
       position='append', action='write', err=900)
else
  open(file=name, unit=nfecra, form='formatted', status='unknown', err=900)
endif

goto 950

900 ierror = 1

950 continue

return
end subroutine
