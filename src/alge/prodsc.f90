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

subroutine prodsc &
!================

 ( ncelet , ncel   , isqrt  , va     , vb     , vavb   )

!===============================================================================
! Purpose:
! --------

! Dot product VAPVB = VA.VB or \/ VA.VB  if ISQRT=1

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! isqrt            ! i  ! <-- ! flag: 1 to return the square root              !
! va, vb(ncelet)   ! ra ! <-- ! vectors to multiply                            !
! vavb             ! r  ! --> ! dot product                                    !
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
double precision vavb
double precision va(ncelet),vb(ncelet)

! Local variables

double precision csdot
external         csdot

!===============================================================================

vavb = csdot(ncel, va, vb)

if (irangp.ge.0) call parsom (vavb)
                 !==========
if (isqrt.eq.1) vavb= sqrt(vavb)

!----
! End
!----

return

end subroutine

