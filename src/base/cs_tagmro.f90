!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
!> \param[in]     dt            time step of the 1D thermal model
!_______________________________________________________________________________

subroutine cs_tagmro &
 ( nfbpcd , ifbpcd ,                          &
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
use elincl
use radiat
use parall
use mesh
use field
use pointe, only:flthr, dflthr
use cs_tagmr

!===============================================================================

implicit none

! Arguments

integer          nfbpcd, ifbpcd(nfbpcd)

double precision dt(ncelet)

! Local variables

integer          ilelt, iel, ifac
integer          ii, jj

double precision dx,dx2
double precision phi,dphi
double precision da(nmur),xa(2,nmur),xsm(nmur)
double precision dtmur(nmur)
double precision tpminf,tpmaxf,tpminp,tpmaxp
double precision dxv, rocp

!===============================================================================

!===============================================================================
! Resolution of the 1-D thermal problem coupled with condensation
!===============================================================================

rocp = rob*cpb

do ilelt = 1, nfbpcd

  ifac = ifbpcd(ilelt)
  iel = ifabor(ifac)

  ! Cote fluide, on recupere le flux et la derivee du flux (implicitation)
  phi = flthr(ilelt)
  dphi= dflthr(ilelt)

  do ii=2, nmur-1

    dx  = dxp(ii)
    dxv = 0.5d0*(dxp(ii-1)+dxp(ii))

    da (ii) = rocp/dt(iel)+theta*condb/(dxp(ii-1)*dxv)            &
                          +theta*condb/(dxp(ii)  *dxv)
    xa(1,ii) =  -theta*condb/(dxp(ii-1)*dxv)
    xa(2,ii) =  -theta*condb/(dxp(ii  )*dxv)
    xsm(ii) = condb*(  tmur(ilelt,ii+1)/(dxp(ii)  *dxv)            &
                     - tmur(ilelt,ii)  /(dxp(ii)  *dxv)            &
                     - tmur(ilelt,ii)  /(dxp(ii-1)*dxv)            &
                     + tmur(ilelt,ii-1)/(dxp(ii-1)*dxv))
  enddo

  ! cote fluide
  ii = 1
  dx  = dxp(ii)
  dx2 = dxp(ii)*dxp(ii)
  da (ii) = rocp/dt(iel)+theta*2.d0*condb/dx2+2.d0*dphi/dx
  xa(1,ii) = 0.d0
  xa(2,ii) = -theta*2.d0*condb/dx2
  xsm(ii) =  2.d0*condb/dx2*(tmur(ilelt,ii+1)-tmur(ilelt,ii))      &
           +(2.d0/dx)*phi

  ! cote externe
  ii = nmur
  dx  = dxp(ii-1)
  dx2 = dxp(ii-1)*dxp(ii-1)
  da (ii) = rocp/dt(iel)+theta*2.d0*condb/dx2+2.d0*hext/dx
  xa(1,ii) = -theta*2.d0*condb/dx2
  xa(2,ii) = 0.d0
  xsm(ii) = 2.d0*condb/dx2*(tmur(ilelt,ii-1)-tmur(ilelt,ii))       &
           -(2.d0/dx)*hext*(tmur(ilelt,ii)-text)

  ! Resolution sur l'increment
  do ii = 1, nmur
    dtmur(ii) = 0.d0
  enddo

  do jj = 1, nmur
    ii = 1
    dtmur(ii) = (xsm(ii)+xa(2,ii)*dtmur(ii+1))/da(ii)
    do ii = 2, nmur-1
      dtmur(ii)= (xsm(ii)+xa(1,ii)*dtmur(ii-1)+xa(2,ii)*dtmur(ii+1)) / da(ii)
    enddo
    ii = nmur
    dtmur(ii) = (xsm(ii)+xa(1,ii)*dtmur(ii-1))/da(ii)
  enddo

  ! Actualisation de la temperature

  do ii = 1, nmur
    tmur(ilelt,ii) = tmur(ilelt,ii)+dtmur(ii)
  enddo

enddo

if (mod(ntcabs,ntlist).eq.0) then

  tpminf = +1.d20
  tpmaxf = -1.d20
  tpminp = +1.d20
  tpmaxp = -1.d20
  do ilelt = 1, nfbpcd
    tpminf = min(tpminf,tmur(ilelt,1))
    tpmaxf = max(tpmaxf,tmur(ilelt,1))
    tpminp = min(tpminp,tmur(ilelt,nmur))
    tpmaxp = max(tpmaxp,tmur(ilelt,nmur))
  enddo

  if (irangp .ge. 0) then
    call parmin(tpminf)
    call parmin(tpminp)
    call parmax(tpmaxf)
    call parmax(tpmaxp)
  endif

  write(nfecra,1000)
  write(nfecra,1001) ttcabs,tpminf,tpmaxf,tpminp,tpmaxp
  write(nfecra,1002)

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
     '------------------------------------'         ,/, &
  3x,' time', 8x,'Tp_f   (min) ',5x,'Tp_f   (max)',6x,  &
                  'Tp_ext(min) ',5x,'Tp_ext (max)'  ,/, &
  3x,'  (s) ',8x, ' (C)       ' ,5x,' (C)        ',6x,  &
                  ' (C)       ' ,5x,' (C)        '  ,/, &
  3x,'------------------------------------------',      &
     '------------------------------------' )
 1001 format( 3x, 5(g15.7,1x) )
 1002 format(&
  3X,'------------------------------------------'   ,   &
     '------------------------------------' )

!----
! End
!----

return

end subroutine cs_tagmro
