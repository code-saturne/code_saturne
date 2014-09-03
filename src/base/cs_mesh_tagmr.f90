!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
! Function:
! ---------

!> \file cs_mesh_tagmr.f90
!>
!> \brief The subroutine is used to generate the 1-D mesh and initialize
!> the temperature field of the thermal model coupled with condensation model.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ifbpcd        index of faces with condensation source terms
!_______________________________________________________________________________

subroutine cs_mesh_tagmr &
 ( nfbpcd , ifbpcd )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use cstphy
use entsor
use ppppar
use ppthch
use ppincl
use elincl
use radiat
use parall
use mesh
use field

use cs_tagmr

!===============================================================================

implicit none

! Arguments

integer          nfbpcd, ifbpcd(nfbpcd)

! Local variables

integer          ilelt, ifac
integer          ii, iter

double precision epsi,delta,r0,r1,epai1

!===============================================================================

!===============================================================================
! 0 - Generate the 1-D mesh of the thermal model coupled
!     with the condensation model in case of solid
!     temperature variable over time
!===============================================================================

if (dxmin .le. 0.d0 .or. dxmin .gt. epais/dble(nmur-1)) then
  !---------------------------------------------------------
  ! Generate a homegeneous 1-D mesh with constant space step
  !---------------------------------------------------------
  do ii = 1,nmur-1
    dxp(ii) = epais/dble(nmur-1)
  enddo
else
  !-----------------------------------------------------------
  ! Generate a heterogeneous 1-D mesh with variable space step
  !-----------------------------------------------------------

  ! Compute the geometric ratio with a iterative method
  iter = 0
  r1   = 2.d0
  delta = 0.0
  epsi = 0.0001d0

  do while (delta .gt. epsi .and. iter .lt. 100)
    iter  = iter + 1
    r0    = r1
    r1    = (1.d0+(epais*(r0-1.d0))/dxmin)**(1.d0/dble(nmur-1))
    epai1 = dxmin*(r1**(nmur-1)-1.d0)/(r1-1.d0)
    delta = abs(epai1-epais)/epais
  end do

  if (iter .ge. 100) then
    write(nfecra,*) '=========================================='
    write(nfecra,*) ' Error with the 1-D mesh Generation       '
    write(nfecra,*) ' Stop inside the cs_mesh_tagmr subroutine '
    write(nfecra,*) '=========================================='
    call csexit(1)
  endif

  ! Comute the final 1-D mesh of the thermal model
  dxp(1)=dxmin
  do ii=2,nmur-1
    dxp(ii) = dxp(ii-1)*r1
  enddo

  write(nfecra,2000) r1
  r0=0.d0
  do ii=1,nmur-1
    r0 = r0 + dxp(ii)
    write(nfecra,2001) ii,dxp(ii),r0
  enddo
  write(nfecra,2002) (dxmin**2)/(2.d0*(condb/(rob*cpb)))
endif

!===============================================================================
! 1 - Initialization of the 1D temperature field for the thermal model
!     which is coupled with the condensation model
!===============================================================================

do ilelt = 1, nfbpcd
  do ii=1, nmur
    ifac = ifbpcd(ilelt)
    tmur(ilelt,ii) = tpar0
  enddo
enddo

!--------
! Formats
!--------

  !================================================================
  2000 format(/,                                                   &
         1x,'=============================================== ',/,  &
         1x,'1-D mesh generation of the thermal model        ',/,  &
            'this one is coupled with the condensation model ',/,  &
         1x,'=============================================== ',/,  &
  4x,' geometric ratio : ',g15.7,/,                                &
  4x,' cell id ',3x,' cell size ',5x,'distance to the wall'   )
 2001 format( 8x,i4,8x,g15.7,6x,g15.7)
 2002 format( /,                                                   &
  4x,'Minimum characteristic time of heat transfer : ',g15.7,/)
  !================================================================

!----
! END
!----

return

end subroutine cs_mesh_tagmr
