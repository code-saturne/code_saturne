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

!> \file yg2xye.f90
!>
!> \brief Compute molar and mass fractions of elementary species Ye, Xe
!>   (fuel, O2, CO2, H2O, N2) from global species Yg (fuel, oxidant, products)
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     yg            global mass fraction
!> \param[out]    ye            elementary mass fraction
!> \param[out]    xe            elementary molar fraction
!_______________________________________________________________________________

subroutine yg2xye (yg, ye, xe)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use cstphy
use entsor
use ppppar
use ppthch
use coincl

!===============================================================================

implicit none

! Arguments

double precision yg(ngazg), ye(ngaze), xe(ngaze)

! Local variables

integer          ige, igg
double precision ytot, mm

!===============================================================================

! Conversion Yg -> Ye
do ige = 1, ngaze
  ye(ige) = 0.d0
  do igg = 1, ngazg
    ye(ige) = ye(ige) + coefeg(ige,igg)*yg(igg)
  enddo
enddo

! Verification
ytot = 0.d0
do ige = 1, ngaze
 ytot = ytot + ye(ige)
enddo

! Warning
if(ytot.lt.0.d0.or.(1.d0-ytot).lt.-epzero) then
  write(nfecra,1000) ytot
endif

! Mixture molar mass molaire
mm = 0.d0
do ige = 1, ngaze
  mm = mm + ye(ige)/wmole(ige)
enddo
mm = 1.d0/mm

! Conversion Ye -> Xe
do ige = 1, ngaze
  xe(ige) = ye(ige)*mm/wmole(ige)
enddo

!--------
! Formats
!--------

1000 format(                                                    /,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : in yg2xye, mass fraction sum exits            ',/,&
'@              physical boundaries [0, 1].                   ',/,&
'@              sum_i=1,ngazge Yi = ',E14.5                    ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return

end subroutine yg2xye
