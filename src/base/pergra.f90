!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine pergra &
!================

 ( ivar   , idimtr , irpvar )

!===============================================================================
! Purpose:
! --------

! Indicate if the variable considered is a component of a vector or tensor
! in the presence of periodicity of rotation

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivar             ! i  ! <-- ! variable number                                !
! idimtr           ! i  ! --> ! 0 if ivar does not match a vector or tensor    !
!                  !    !     !   or there is no periodicity of rotation       !
!                  !    !     ! 1 for velocity, 2 for Reynolds stress          !
! irpvar           ! i  ! --> ! -1 if ivar does not match a vector or tensor   !
!                  !    !     ! In presence of periodicity of rotation:        !
!                  !    !     !  0 for iu, 1 for iv, 2 for iw                  !
!                  !    !     !  0 for ir11, 1 for ir22, 2 for ir33,           !
!                  !    !     !  3 for ir12, 4 for ir13, 5 for ir23            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use period

!===============================================================================

implicit none

! Arguments

integer, intent(in)  :: ivar
integer, intent(out) :: idimtr, irpvar

! Local variables

!===============================================================================

idimtr = 0
irpvar = -1

if (iperot .eq. 0) return

if (ivar.eq. iu) then
  idimtr = 1
  irpvar = 0
else if (ivar.eq. iv) then
  idimtr = 1
  irpvar = 1
else if (ivar.eq. iw) then
  idimtr = 1
  irpvar = 2
else if (itytur.eq.3) then
  if (ivar.eq. ir11) then
    irpvar = 0
  else if (ivar.eq. ir22) then
    irpvar = 1
  else if (ivar.eq. ir33) then
    irpvar = 2
  else if (ivar.eq. ir12) then
    irpvar = 3
  else if (ivar.eq. ir13) then
    irpvar = 4
  else if (ivar.eq. ir23) then
    irpvar = 5
  endif
  if (irpvar .ge. 0) idimtr = 2
endif

!----
! End
!----

return
end subroutine
