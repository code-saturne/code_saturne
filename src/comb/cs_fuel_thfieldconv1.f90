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
! --------
!> \file cs_fuel_thfieldconv1.f90
!> \brief Calculation of the gas temperature
!>  Function with the gas enthalpy and concentrations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     eh            gas enthalpy
!>                              (\f$ j . kg \f$ of gaseous mixture)
!> \param[in]     fuel1         mass fraction CHx1
!> \param[in]     fuel2         mass fraction CHx2
!> \param[in]     fuel3         mass fraction CO
!> \param[in]     oxyd          mass fraction O2
!> \param[in]     prod1         mass fraction CO2
!> \param[in]     prod2         mass fraction H2O
!> \param[in]     xiner         mass fraction N2
!> \param[in,out] tp            gas temperature (in kelvin)
!______________________________________________________________________________!

subroutine cs_fuel_thfieldconv1 &
 ( ncelet , ncel   ,                                              &
   eh     ,                                                       &
   fuel1  , fuel2  , fuel3  , fuel4 , fuel5 , fuel6 , fuel7 ,     &
   oxyd   , prod1  , prod2  , prod3 , xiner ,                     &
   tp     )

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision eh(ncelet)
double precision fuel1(ncelet), fuel2(ncelet) , fuel3(ncelet)
double precision fuel4(ncelet), fuel5(ncelet) , fuel6(ncelet)
double precision fuel7(ncelet)
double precision oxyd(ncelet), prod1(ncelet),prod2(ncelet)
double precision prod3(ncelet),xiner(ncelet)
double precision tp(ncelet)

! Local variables

integer          ii, icel
double precision eh0,eh1

!===============================================================================

ii = npo-1
do icel = 1, ncel

  ! --- Eventual clipping of temperature at TH(NPO) if EH > EH1

  eh1 = fuel1(icel)*ehgaze(ifo0,ii+1)                       &
       +fuel2(icel)*ehgaze(ifov,ii+1)                       &
       +fuel3(icel)*ehgaze(ico ,ii+1)                       &
       +fuel4(icel)*ehgaze(ih2s,ii+1)                       &
       +fuel5(icel)*ehgaze(ihy ,ii+1)                       &
       +fuel6(icel)*ehgaze(ihcn,ii+1)                       &
       +fuel7(icel)*ehgaze(inh3,ii+1)                       &
       +oxyd(icel) *ehgaze(io2 ,ii+1)                       &
       +prod1(icel)*ehgaze(ico2,ii+1)                       &
       +prod2(icel)*ehgaze(ih2o,ii+1)                       &
       +prod3(icel)*ehgaze(iso2,ii+1)                       &
       +xiner(icel)*ehgaze(in2 ,ii+1)
  if ( eh(icel) .ge. eh1 ) tp(icel)= th(ii+1)
enddo

ii = 1
do icel = 1, ncel

  ! --- Eventual clipping of temperature at TH(1) if EH < EH0

  eh0= fuel1(icel)*ehgaze(ifo0,ii)                         &
      +fuel2(icel)*ehgaze(ifov,ii)                         &
      +fuel3(icel)*ehgaze(ico ,ii)                         &
      +fuel4(icel)*ehgaze(ih2s,ii)                         &
      +fuel5(icel)*ehgaze(ihy ,ii)                         &
      +fuel6(icel)*ehgaze(ihcn,ii)                         &
      +fuel7(icel)*ehgaze(inh3,ii)                         &
      +oxyd(icel) *ehgaze(io2 ,ii)                         &
      +prod1(icel)*ehgaze(ico2,ii)                         &
      +prod2(icel)*ehgaze(ih2o,ii)                         &
      +prod3(icel)*ehgaze(iso2,ii)                         &
      +xiner(icel)*ehgaze(in2 ,ii)

  if ( eh(icel) .le. eh0 ) then
    tp(icel)= th(1)
  endif
enddo

do ii = 1, npo-1
  do icel = 1, ncel

  eh0 = fuel1(icel)*ehgaze(ifo0,ii)                        &
       +fuel2(icel)*ehgaze(ifov,ii)                        &
       +fuel3(icel)*ehgaze(ico ,ii)                        &
       +fuel4(icel)*ehgaze(ih2s,ii)                        &
       +fuel5(icel)*ehgaze(ihy ,ii)                        &
       +fuel6(icel)*ehgaze(ihcn,ii)                        &
       +fuel7(icel)*ehgaze(inh3,ii)                        &
       +oxyd(icel) *ehgaze(io2 ,ii)                        &
       +prod1(icel)*ehgaze(ico2,ii)                        &
       +prod2(icel)*ehgaze(ih2o,ii)                        &
       +prod3(icel)*ehgaze(iso2,ii)                        &
       +xiner(icel)*ehgaze(in2 ,ii)

  eh1 = fuel1(icel)*ehgaze(ifo0,ii+1)                      &
       +fuel2(icel)*ehgaze(ifov,ii+1)                      &
       +fuel3(icel)*ehgaze(ico ,ii+1)                      &
       +fuel4(icel)*ehgaze(ih2s,ii+1)                      &
       +fuel5(icel)*ehgaze(ihy ,ii+1)                      &
       +fuel6(icel)*ehgaze(ihcn,ii+1)                      &
       +fuel7(icel)*ehgaze(inh3,ii+1)                      &
       +oxyd(icel) *ehgaze(io2 ,ii+1)                      &
       +prod1(icel)*ehgaze(ico2,ii+1)                      &
       +prod2(icel)*ehgaze(ih2o,ii+1)                      &
       +prod3(icel)*ehgaze(iso2,ii+1)                      &
       +xiner(icel)*ehgaze(in2 ,ii+1)

    if ( eh(icel) .ge. eh0 .and. eh(icel) .le. eh1  ) then
      tp(icel)= th(ii) + (eh(icel)-eh0)                    &
                       *(th(ii+1)-th(ii))/(eh1-eh0)
    endif
  enddo
enddo

!----
! End
!----

return
end subroutine
