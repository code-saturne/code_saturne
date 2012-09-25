!-------------------------------------------------------------------------------

!VERS

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

subroutine usporo
!================

!===============================================================================
! Function :
! ----------
! Compute the porosity (volume factor) when module is activated (iporos = 1)
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use pointe
use cstnum
use parall
use period
use mesh

!===============================================================================

implicit none

! Local variables

integer          iel, ii, jj
double precision x, pormin, pormax, hc, ll, dhc


!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
!===============================================================================

! Example: fixe a linear porosity profil

do iel = 1, ncel
  x = xyzcen(1,iel)
  if (x.le.(ll/2.d0)) then
    hc = 1.d0 - 2.d0*dhc*x/ll
  else
    hc = 1.d0 - dhc + dhc*(2.d0*x-ll)/ll
  endif

  porosi(iel) = hc

! TODO move elsewhere
!  if (porosi(iel).lt.0.d0) then
!    write(nfecra,*) 'Negative porosity'
!    call csexit(1)
!  elseif (porosi(iel).gt.1.d0) then
!    write(nfecra,*) 'Porosity stricly greater than 1.0'
!    call csexit(1)
!  endif

  pormin = min(pormin,porosi(iel))
  pormax = max(pormax,porosi(iel))
enddo

! Periodicity and parallelism treatment
if(irangp.ge.0) then
  call parmax (pormax)
  call parmin (pormin)
endif

if (iperio.eq.1.or.irangp.ge.0) then
  call synsca(porosi)
endif

return
end subroutine usporo
