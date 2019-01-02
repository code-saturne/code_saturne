!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file d3phst.f90
!>
!> \brief Specific physic subroutine: diffusion flame.
!>
!> Calculation of local stoechiometric enthalpy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncel          number of cells
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     indpdf        indicator for pdf integration or mean value
!> \param[in]     dirmin        Dirac's peak value at \f$ f_{min} \f$
!> \param[in]     dirmax        Dirac's peak value at \f$ f_{max} \f$
!> \param[in]     fdeb          abscissa of rectangle low boundary
!> \param[in]     ffin          abscissa of rectangle high boundary
!> \param[in]     hrec          rectangle height
!> \param[in]     fm            mean mixture fraction at cell centers
!> \param[in]     hm            mean mixture enthalpy at cell centers
!> \param[in,out] hstoe         local stoechiometric enthalpy
!_______________________________________________________________________________

subroutine d3phst &
 ( ncelet , ncel   , indpdf ,                                     &
   dirmin , dirmax , fdeb   , ffin   , hrec   ,                   &
   fm     , hm     ,                                              &
   hstoe  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          indpdf(ncelet)
double precision dirmin(ncelet), dirmax(ncelet)
double precision fdeb(ncelet), ffin(ncelet), hrec(ncelet)
double precision fm(ncelet), hm(ncelet), hstoe(ncelet)


! Local variables

integer          iel
double precision fsir, hhh, hct, f1, f2
double precision epsi

integer          n1, n2
double precision hsmax, hsmin


!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

epsi  = 1.d-6
fsir =  fs(1)

n1 = 0
n2 = 0
hsmin = grand
hsmax =-grand


do iel = 1, ncel

  if ( indpdf(iel) .eq. 0 ) then

!===============================================================================
! 1. DETERMINATION DE HSTOE SANS INTEGRATION
!===============================================================================

    if (fm(iel).le.fsir .and. fm(iel).gt.epsi) then
      hstoe(iel) = (fsir*hm(iel)+(fm(iel)-fsir)*hinoxy)        &
                  / fm(iel)
    elseif (fm(iel).lt.(1.d0-epsi)) then
      hstoe(iel) = ((1.d0-fsir)*hm(iel)+(fsir-fm(iel))*hinfue) &
                 / (1.d0-fm(iel))
    endif

  else

!===============================================================================
! 2. DETERMINATION DE HSTOE AVEC INTEGRATION
!===============================================================================

    hct = dirmin(iel)*hinoxy+dirmax(iel)*hinfue
    hhh = 0.d0
    if (hrec(iel).gt.epsi) then

      if (fdeb(iel).le.fsir) then
        f1 = fdeb(iel)
        f2 = min(fsir,ffin(iel))
        hct = hct + hrec(iel)*                                   &
              (f2-f1)*hinoxy*(2.d0*fsir-f1-f2)/(2.d0*fsir)
        hhh = hhh + hrec(iel)*(f2**2-f1**2)/(2.d0*fsir)
      endif
      if (ffin(iel).gt.fsir) then
        f1 = max( fsir,fdeb(iel))
        f2 = ffin(iel)
        hct = hct + hrec(iel) *                                  &
             (f2-f1)*hinfue*(f2+f1-2.d0*fsir)/(2.d0*(1.d0-fsir))
        hhh = hhh +                                               &
              hrec(iel)*(f2-f1)*(2.d0-f1-f2)/(2.d0*(1.d0-fsir))
      endif
      hstoe(iel) = (hm(iel)-hct)/ hhh

    endif

  endif

  ! Clipping a HSTOEA = HH(1)     en max
  ! Clipping a HSTOEA = HH(NMAXH) em min

  if (hstoe(iel) .gt. hh(1)) then
    n1 = n1 + 1
    hsmax = max(hstoe(iel),hsmax)
    hstoe(iel) = hh(1)
  endif

  if (hstoe(iel) .lt. hh(nmaxh)) then
    n2 = n2 + 1
    hsmin = min(hstoe(iel),hsmin)
    hstoe(iel) = hh(nmaxh)
  endif

enddo

if (irangp.ge.0) then
  call parcpt (n1)
  !==========
  call parcpt (n2)
  !==========
  call parmax (hsmax)
  !==========
  call parmin (hsmin)
  !==========
endif

if ( n1.gt.0 ) then
  write(nfecra,1000) n1,hsmax,hh(1)
endif
if ( n2.gt.0 ) then
  write(nfecra,1001) n2,hsmin,hh(nmaxh)
endif

!----
! FORMATS
!----

 1000   format(1X,' Clipping de HSTOE EN MAX EN ',I8,' POINTS',/, &
         1X,'     Valeur Max : ',G15.7,/,                   &
         1X,'     Valeur De Clipping : ',G15.7,/)
 1001   format(1X,' Clipping de HSTOE EN MIN EN ',I8,' POINTS',/, &
         1X,'     Valeur Max : ',G15.7,/,                   &
         1X,'     Valeur De Clipping : ',G15.7,/)

return
end subroutine

