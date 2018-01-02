!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
! Function :
! --------

!> \file distpr2.f90
!> \brief Compute distance to wall by a brute force geometric approach
!>        (serial only)
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     itypfb        boundary face types
!> \param[out]    distpa        distance to wall
!______________________________________________________________________________

subroutine distpr2 &
 ( itypfb , distpa )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use mesh
use parall
use period
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)

double precision distpa(ncelet)

! Local variables

integer          ifac  , iel

double precision xdis, dismax, dismin

!===============================================================================

! normalement, on ne passe pas en parallele ici,  mais au cas ou ...
if (irangp.ge.0 .or. iperio.gt.0) then
  call csexit(1)
endif

!===============================================================================
! Deprecated model to compute wall distance
!===============================================================================

! on fera attention en parallelisme ou periodicite
!    (une paroi peut etre plus proche en traversant un bord ...)

do iel = 1, ncel
  distpa(iel) = grand*grand
enddo

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    do iel = 1, ncel
      xdis =   (cdgfbo(1,ifac)-xyzcen(1,iel))**2            &
             + (cdgfbo(2,ifac)-xyzcen(2,iel))**2            &
             + (cdgfbo(3,ifac)-xyzcen(3,iel))**2
      if (distpa(iel).gt.xdis) then
        distpa(iel) = xdis
      endif
    enddo
  endif
enddo

do iel = 1, ncel
  distpa(iel) = sqrt(distpa(iel))
enddo

!===============================================================================
! Compute bounds and print info
!===============================================================================

dismax = -grand
dismin =  grand

do iel = 1, ncel
  dismin = min(distpa(iel),dismin)
  dismax = max(distpa(iel),dismax)
enddo

write(nfecra,1000) dismin, dismax

!===============================================================================
! Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI                                      ',/,&
'    -------------------                                      ',/,&
'                                                             ',/,&
'   Distance min = ',E14.5    ,'  Distance max = ',E14.5      ,/)

#else

 1000 format(                                                           &
'                                                             ',/,&
' ** WALL DISTANCE                                            ',/,&
'    -------------                                            ',/,&
'                                                             ',/,&
'  Min distance = ',E14.5    ,' Max distance = ',E14.5         ,/)

#endif

!----
! End
!----

return
end subroutine distpr2
