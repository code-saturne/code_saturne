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

!===============================================================================
! Function:
! ---------

!> \file cs_tagmro.f90
!>
!> \brief The 1D thermal model to compute the temperature to impose
!> at the cold wall. This one is used by the COPAIN model to estimate
!> the heat flux at the wall where the condensation occurs.
!>
!> This subroutine is used to compute at each face the
!> \f$T^{fb}_{\mbox{mur}} \f$ at cold wall.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     izzftcd        faces zone with condensation source terms imposed
!>                              (at previous and current time steps)
!> \param[in]     dt            time step of the 1D thermal model
!_______________________________________________________________________________

subroutine cs_tagmro &
 ( nfbpcd , ifbpcd , izzftcd ,                    &
   dt     )

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
use pointe, only:flthr, dflthr
use cs_nz_condensation, only:nzones,iztag1d
use cs_nz_tagmr

!===============================================================================

implicit none

! Arguments

integer          nfbpcd, ifbpcd(nfbpcd), izzftcd(nfbpcd)

double precision dt(ncelet)

! Local variables

integer          ii, iel, ifac
integer          kk, iz

double precision dx,dx2
double precision phi,dphi
double precision da(znmurx),xa(2,znmurx),xsm(znmurx)
double precision dtmur(znmurx)
double precision tpminf(nzones),tpmaxf(nzones),tpminp(nzones),tpmaxp(nzones)
double precision dxv, rocp

!===============================================================================

!===============================================================================
! Resolution of the 1-D thermal problem coupled with condensation
!===============================================================================


do ii = 1, nfbpcd

  ifac = ifbpcd(ii)
  iel  = ifabor(ifac)

  iz = izzftcd(ii)

  if(iztag1d(iz).eq.1) then

    rocp = zrob(iz)*zcpb(iz)
    ! Cote fluide, on recupere le flux et la derivee du flux (implicitation)
    phi = flthr(ii)
    dphi= dflthr(ii)

    do kk=2, znmur(iz)-1

      dx  = zdxp(iz,kk)
      dxv = 0.5d0*(zdxp(iz,kk-1)+zdxp(iz,kk))

      da (kk) = rocp/dt(iel)+ztheta(iz)*zcondb(iz)/(zdxp(iz,kk-1)*dxv)      &
                            +ztheta(iz)*zcondb(iz)/(zdxp(iz,kk)  *dxv)
      xa(1,kk) =  -ztheta(iz)*zcondb(iz)/(zdxp(iz,kk-1)*dxv)
      xa(2,kk) =  -ztheta(iz)*zcondb(iz)/(zdxp(iz,kk  )*dxv)
      xsm(kk) = zcondb(iz)*( ztmur(ii,kk+1)/(zdxp(iz,kk)  *dxv)            &
                       - ztmur(ii,kk)  /(zdxp(iz,kk)  *dxv)                &
                       - ztmur(ii,kk)  /(zdxp(iz,kk-1)*dxv)                &
                       + ztmur(ii,kk-1)/(zdxp(iz,kk-1)*dxv))
    enddo

    ! cote fluide
    kk = 1
    dx  = zdxp(iz,kk)
    dx2 = zdxp(iz,kk)*zdxp(iz,kk)
    da (kk)  = rocp/dt(iel)+ztheta(iz)*2.d0*zcondb(iz)/dx2+2.d0*dphi/dx
    xa(1,kk) = 0.d0
    xa(2,kk) = -ztheta(iz)*2.d0*zcondb(iz)/dx2
    xsm(kk)  =  2.d0*zcondb(iz)/dx2*(ztmur(ii,kk+1)-ztmur(ii,kk))      &
             +(2.d0/dx)*phi

    ! cote externe
    kk = znmur(iz)
    dx  = zdxp(iz,kk-1)
    dx2 = zdxp(iz,kk-1)*zdxp(iz,kk-1)
    da (kk)  = rocp/dt(iel)+ztheta(iz)*2.d0*zcondb(iz)/dx2+2.d0*zhext(iz)/dx
    xa(1,kk) = -ztheta(iz)*2.d0*zcondb(iz)/dx2
    xa(2,kk) = 0.d0
    xsm(kk)  = 2.d0*zcondb(iz)/dx2*(ztmur(ii,kk-1)-ztmur(ii,kk))       &
             -(2.d0/dx)*zhext(iz)*(ztmur(ii,kk)-ztext(iz))

    ! Resolution sur l'increment
    do kk = 1, znmur(iz)
      dtmur(kk) = 0.d0
    enddo

    kk = 1
    dtmur(kk) = (xsm(kk)+xa(2,kk)*dtmur(kk+1))/da(kk)
    do kk = 2, znmur(iz)-1
      dtmur(kk)= (xsm(kk)+xa(1,kk)*dtmur(kk-1)+xa(2,kk)*dtmur(kk+1)) / da(kk)
    enddo
    kk = znmur(iz)
    dtmur(kk) = (xsm(kk)+xa(1,kk)*dtmur(kk-1))/da(kk)

    ! Actualisation de la temperature

    do kk = 1, znmur(iz)
      ztmur(ii,kk) = ztmur(ii,kk)+dtmur(kk)
    enddo

  endif
enddo

if (mod(ntcabs,ntlist).eq.0) then

  do iz = 1, nzones
    tpminf(iz) = +1.d20
    tpmaxf(iz) = -1.d20
    tpminp(iz) = +1.d20
    tpmaxp(iz) = -1.d20
  enddo

  do ii = 1, nfbpcd
    iz = izzftcd(ii)
    if(iztag1d(iz).eq.1) then
      tpminf(iz) = min(tpminf(iz),ztmur(ii,1))
      tpmaxf(iz) = max(tpmaxf(iz),ztmur(ii,1))
      tpminp(iz) = min(tpminp(iz),ztmur(ii,znmur(iz)))
      tpmaxp(iz) = max(tpmaxp(iz),ztmur(ii,znmur(iz)))
    endif
  enddo

  do iz = 1, nzones
    if (irangp.ge.0) then
      call parmin(tpminf(iz))
      call parmin(tpminp(iz))
      call parmax(tpmaxf(iz))
      call parmax(tpmaxp(iz))
    endif

    if (irangp.le.0) then
      write(nfecra,1000)
      write(nfecra,1001) ttcabs,iz,tpminf(iz),tpmaxf(iz),tpminp(iz),tpmaxp(iz)
      write(nfecra,1002)
    endif
  enddo

endif

!--------
! Formats
!--------

  1000 format(/,&
         3x,'===================================== ',/, &
         3x,'Resolution of the 1-D thermal problem ',/, &
         3x,' coupled with the condensation model  ',/, &
         3x,'===================================== ',/, &
              /,&
  3x,'------------------------------------------'   ,   &
     '----------------------------------------------',/,&
  3x,' time', 8x, ' izones', 5x,&
                  'Tp_f  (min) ',5x,'Tp_f   (max)',6x,  &
                  'Tp_ext(min) ',5x,'Tp_ext (max)'  ,/, &
  3x,'  (s) ',8x, '       ', 5x,&
                  ' (C)       ' ,5x,' (C)        ',6x,  &
                  ' (C)       ' ,5x,' (C)        '  ,/, &
  3x,'------------------------------------------',      &
     '----------------------------------------------' )
 1001 format( 3x, g15.7, 1x, i4, 3x, 4(g15.7,1x) )
 1002 format(&
  3X,'------------------------------------------'   ,   &
     '----------------------------------------------' )

!----
! End
!----

return

end subroutine cs_tagmro
