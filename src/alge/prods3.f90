!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine prods3 &
!================

 ( ncelet , ncel   , isqrt  ,                                     &
   va1    , vb1    , va2    , vb2    , va3    , vb3    ,          &
   vavb1  , vavb2  , vavb3  )

!===============================================================================
! Purpose:
! --------

! Triple dot product VAPVB = VA.VB or \/ VA.VB  if ISQRT=1

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! isqrt            ! i  ! <-- ! flag: 1 to return the square root              !
! va1(), vb1()     ! ra ! <-- ! first vectors to multiply                      !
! va2(), vb2()     ! ra ! <-- ! second vectors to multiply                     !
! va3(), vb3()     ! ra ! <-- ! third vectors to multiply                      !
! vavb1            ! r  ! --> ! first dot product                              !
! vavb2            ! r  ! --> ! second dot product                             !
! vavb3            ! r  ! --> ! third dot product                              !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet,ncel,isqrt
double precision vavb1, vavb2, vavb3
double precision va1(ncelet),vb1(ncelet)
double precision va2(ncelet),vb2(ncelet)
double precision va3(ncelet),vb3(ncelet)

! Local variables

integer nvavb
double precision vavb(3)

double precision csdot
external         csdot

!===============================================================================

vavb(1) = csdot(ncel, va1, vb1)
vavb(2) = csdot(ncel, va2, vb2)
vavb(3) = csdot(ncel, va3, vb3)

if (irangp.ge.0) then
  nvavb = 3
  call parrsm (nvavb, vavb)
  !==========
endif

vavb1 = vavb(1)
vavb2 = vavb(2)
vavb3 = vavb(3)

if (isqrt.eq.1) then
  vavb1 = sqrt(vavb1)
  vavb2 = sqrt(vavb2)
  vavb3 = sqrt(vavb3)
endif

!----
! End
!----

return

end subroutine

