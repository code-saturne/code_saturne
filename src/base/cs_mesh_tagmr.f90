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
!> \param[in]     izzftcd       faces zone with condensation source terms imposed
!>                              (at previous and current time steps)
!_______________________________________________________________________________

subroutine cs_mesh_tagmr &
 ( nfbpcd , izzftcd )

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
use radiat
use parall
use mesh
use field

use cs_nz_tagmr

!===============================================================================

implicit none

! Arguments

integer          nfbpcd, izzftcd(nfbpcd)

! Local variables

integer          ii, kk, iz
integer          iter

double precision epsi,delta,r0,r1,epai1

!===============================================================================

!===============================================================================
! 0 - Generate the 1-D mesh of the thermal model coupled
!     with the condensation model in case of solid
!     temperature variable over time
!===============================================================================

do ii = 1, nfbpcd

  iz  = izzftcd(ii)
  if (zdxmin(iz).le.0.d0 .or. zdxmin(iz).gt.zepais(iz)/dble(znmur(iz)-1)) then
    !---------------------------------------------------------
    ! Generate a homegeneous 1-D mesh with constant space step
    !---------------------------------------------------------
    do kk = 1, znmur(iz)-1
      zdxp(iz,kk) = zepais(iz)/dble(znmur(iz)-1)
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
      r1    = (1.d0+(zepais(iz)*(r0-1.d0))/zdxmin(iz))**(1.d0/dble(znmur(iz)-1))
      epai1 = zdxmin(iz)*(r1**(znmur(iz)-1)-1.d0)/(r1-1.d0)
      delta = abs(epai1-zepais(iz))/zepais(iz)
    end do

    if (iter .ge. 100) then
      write(nfecra,*) '=========================================='
      write(nfecra,*) ' Error with the 1-D mesh Generation       '
      write(nfecra,*) ' Stop inside the cs_mesh_tagmr subroutine '
      write(nfecra,*) '=========================================='
      call csexit(1)
    endif

    ! Comute the final 1-D mesh of the thermal model
    zdxp(iz,1)= zdxmin(iz)
    do kk = 2, znmur(iz)-1
      zdxp(iz,kk) = zdxp(iz,kk-1)*r1
    enddo

    write(nfecra,2000) r1
    r0=0.d0
    do kk = 1, znmur(iz)-1
      r0 = r0 + zdxp(iz,kk)
      write(nfecra,2001) kk,zdxp(iz,kk),r0
    enddo
    write(nfecra,2002) (zdxmin(iz)**2)/(2.d0*(zcondb(iz)/(zrob(iz)*zcpb(iz))))
  endif

enddo
!===============================================================================
! 1 - Initialization of the 1D temperature field for the thermal model
!     which is coupled with the condensation model
!===============================================================================

if(isuite.eq.0) then
  do ii = 1, nfbpcd
    iz = izzftcd(ii)
    do kk = 1, znmur(iz)
      ztmur(ii,kk) = ztpar0(iz)
    enddo
  enddo
endif
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
