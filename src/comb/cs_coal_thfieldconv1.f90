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
!> \file cs_coal thfieldconv1.f90
!> \brief Calculation of the gas temperature
!>        Function with gas enthalpy and concentrations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet          number of extended (real + ghost) cells
!> \param[in]     ncel            number of cells
!> \param[in]     ntbmci          macro table size mc integers
!> \param[in]     ntbmcr          macro table size mc reals
!> \param[in]     eh              gas enthalpy
!>                                (j/kg of gaseous mixture)
!> \param[in]     fuel1           mass fraction CHx1
!> \param[in]     fuel2           mass fraction CHx2
!> \param[in]     fuel3           mass fraction CO
!> \param[in]     oxyd            mass fraction O2
!> \param[in]     prod1           mass fraction CO2
!> \param[in]     prod2           mass fraction H2O
!> \param[in]     xiner           mass fraction N2
!> \param[in,out] tp              gas temperature in kelvin
!> \param[in]     tbmci           macro integer table mc travail
!> \param[in]     tbmcr           macro real table    mc travail
!> \param[in,out] eh0             real work array
!> \param[in,out] eh1             real work array
!______________________________________________________________________________!

subroutine cs_coal_thfieldconv1 &
 ( ncelet , ncel   ,                                              &
   eh     ,          rtp    ,                                     &
   fuel1  , fuel2  , fuel3  , fuel4 , fuel5 , fuel6 , fuel7 ,     &
   oxyd   , prod1  , prod2  , prod3 , xiner ,                     &
   tp     )

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use dimens, only: nvar
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision eh(ncelet)   , rtp(ncelet,nflown:nvar)
double precision fuel1(ncelet), fuel2(ncelet), fuel3(ncelet)
double precision fuel4(ncelet), fuel5(ncelet), fuel6(ncelet) , fuel7(ncelet)
double precision oxyd(ncelet), xiner(ncelet)
double precision prod1(ncelet),prod2(ncelet),prod3(ncelet)
double precision tp(ncelet)

! Local variables

integer          i, icel, icha

double precision ychx10 , ychx20 , ehchx1 , ehchx2
double precision den1   , den2 , eh0 , eh1

integer          iok
double precision , dimension ( : , : )     , allocatable :: f1mc,f2mc
double precision, dimension(:), pointer :: x1

!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", x1)

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(f1mc(1:ncel,1:ncharb),f2mc(1:ncel,1:ncharb),stat=iok)
!----
if ( iok > 0  ) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '    cs_coal_thfieldconv1         '
  call csexit(1)
endif
!===============================================================================

do icel = 1, ncel
  do icha = 1, ncharb
    f1mc(icel,icha) = rtp(icel,isca(if1m(icha))) / x1(icel)
    f2mc(icel,icha) = rtp(icel,isca(if2m(icha))) / x1(icel)
  enddo
enddo

i = npo-1
do icel = 1, ncel

! --- Calculation of enthalpy of the gaseous species CHx1m
!                                            and CHx2m at TH(NPO)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb

    den1   = 1.d0                                                  &
         / ( a1(icha)*wmole(ichx1c(icha))                          &
            +b1(icha)*wmole(ico)                                   &
            +c1(icha)*wmole(ih2o)                                  &
            +d1(icha)*wmole(ih2s)                                  &
            +e1(icha)*wmole(ihcn)                                  &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10                                                &
            +den1*(f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1                                                &
            +den1*( ehgaze(ichx1c(icha),i+1)                       &
                   *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                  &
         / ( a2(icha)*wmole(ichx2c(icha))                          &
            +b2(icha)*wmole(ico)                                   &
            +c2(icha)*wmole(ih2o)                                  &
            +d2(icha)*wmole(ih2s)                                  &
            +e2(icha)*wmole(ihcn)                                  &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20                                                &
            +den2 *(f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                       &
         ( ehgaze(ichx2c(icha),i+1)                                &
          *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i+1)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i+1)
  endif

  ! --- Eventual clipping of temperature at TH(NPO) if EH > EH1

  eh1 = fuel1(icel)*ehchx1                                  &
       +fuel2(icel)*ehchx2                                  &
       +fuel3(icel)*ehgaze(ico ,i+1)                        &
       +fuel4(icel)*ehgaze(ih2s,i+1)                        &
       +fuel5(icel)*ehgaze(ihy ,i+1)                        &
       +fuel6(icel)*ehgaze(ihcn,i+1)                        &
       +fuel7(icel)*ehgaze(inh3,i+1)                        &
       +oxyd(icel) *ehgaze(io2 ,i+1)                        &
       +prod1(icel)*ehgaze(ico2,i+1)                        &
       +prod2(icel)*ehgaze(ih2o,i+1)                        &
       +prod3(icel)*ehgaze(iso2,i+1)                        &
       +xiner(icel)*ehgaze(in2 ,i+1)

  if ( eh(icel).ge.eh1 ) tp(icel)= th(i+1)
enddo

i = 1
do icel = 1, ncel

  ! --- Calculation of enthalpy of the gaseous species CHx1m
  !                                            and CHx2m at TH(1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                    &
         / ( a1(icha)*wmole(ichx1c(icha))                            &
            +b1(icha)*wmole(ico)                                     &
            +c1(icha)*wmole(ih2o)                                    &
            +d1(icha)*wmole(ih2s)                                    &
            +e1(icha)*wmole(ihcn)                                    &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10                                                  &
            +den1*(f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1                                                  &
            +den1*( ehgaze(ichx1c(icha),i)                           &
                   *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                     &
         / ( a2(icha)*wmole(ichx2c(icha))                             &
            +b2(icha)*wmole(ico)                                      &
            +c2(icha)*wmole(ih2o)                                     &
            +d2(icha)*wmole(ih2s)                                     &
            +e2(icha)*wmole(ihcn)                                     &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20                                                   &
            +den2*(f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2                                                   &
            +den2*( ehgaze(ichx2c(icha),i)                            &
                   *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i)
  endif

  ! --- Eventual clipping of temperature at TH(1) if EH < EH0

  eh0  = fuel1(icel)*ehchx1                                  &
        +fuel2(icel)*ehchx2                                  &
        +fuel3(icel)*ehgaze(ico ,i)                          &
        +fuel4(icel)*ehgaze(ih2s,i)                          &
        +fuel5(icel)*ehgaze(ihy ,i)                          &
        +fuel6(icel)*ehgaze(ihcn,i)                          &
        +fuel7(icel)*ehgaze(inh3,i)                          &
        +oxyd(icel) *ehgaze(io2 ,i)                          &
        +prod1(icel)*ehgaze(ico2,i)                          &
        +prod2(icel)*ehgaze(ih2o,i)                          &
        +prod3(icel)*ehgaze(iso2,i)                          &
        +xiner(icel)*ehgaze(in2 ,i)

  if ( eh(icel) .le. eh0 ) tp(icel)= th(i)

enddo


do i = 1, npo-1
  do icel = 1, ncel

    ! --- Calculation of enthalpy of the gaseous species CHx1m
    !                                            and CHx2m for TH(I)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                                 &
             / ( a1(icha)*wmole(ichx1c(icha))                       &
                +b1(icha)*wmole(ico)                                &
                +c1(icha)*wmole(ih2o)                               &
                +d1(icha)*wmole(ih2s)                               &
                +e1(icha)*wmole(ihcn)                               &
                +f1(icha)*wmole(inh3) )
      ychx10 = ychx10                                               &
              +den1*f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha))
      ehchx1 = ehchx1                                               &
              +den1*( ehgaze(ichx1c(icha),i)                        &
                     *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                                 &
             / ( a2(icha)*wmole(ichx2c(icha))                       &
                +b2(icha)*wmole(ico)                                &
                +c2(icha)*wmole(ih2o)                               &
                +d2(icha)*wmole(ih2s)                               &
                +e2(icha)*wmole(ihcn)                               &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20                                               &
              +den2*(f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2                                               &
              +den2*( ehgaze(ichx2c(icha),i)                        &
                     *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i)
    endif
    eh0 = fuel1(icel)*ehchx1                                  &
         +fuel2(icel)*ehchx2                                  &
         +fuel3(icel)*ehgaze(ico ,i)                          &
         +fuel4(icel)*ehgaze(ih2s,i)                          &
         +fuel5(icel)*ehgaze(ihy ,i)                          &
         +fuel6(icel)*ehgaze(ihcn,i)                          &
         +fuel7(icel)*ehgaze(inh3,i)                          &
         +oxyd(icel) *ehgaze(io2 ,i)                          &
         +prod1(icel)*ehgaze(ico2,i)                          &
         +prod2(icel)*ehgaze(ih2o,i)                          &
         +prod3(icel)*ehgaze(iso2,i)                          &
         +xiner(icel)*ehgaze(in2 ,i)

    ! --- Calculation of enthalpy of the gaseous species CHx1m
    !                                            and CHx2m for TH(I+1)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                                   &
           / ( a1(icha)*wmole(ichx1c(icha))                           &
              +b1(icha)*wmole(ico)                                    &
              +c1(icha)*wmole(ih2o)                                   &
              +d1(icha)*wmole(ih2s)                                   &
              +e1(icha)*wmole(ihcn)                                   &
              +f1(icha)*wmole(inh3) )
      ychx10 = ychx10                                                 &
              +den1*f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha))
      ehchx1 = ehchx1                                                 &
              +den1*( ehgaze(ichx1c(icha),i+1)                        &
                     *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                                   &
             / ( a2(icha)*wmole(ichx2c(icha))                         &
                +b2(icha)*wmole(ico)                                  &
                +c2(icha)*wmole(ih2o)                                 &
                +d2(icha)*wmole(ih2s)                                 &
                +e2(icha)*wmole(ihcn)                                 &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20                                                 &
              +den2*f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha))
      ehchx2 = ehchx2                                                 &
              +den2*( ehgaze(ichx2c(icha),i+1)                        &
                     *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i+1)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i+1)
    endif

    eh1 = fuel1(icel)*ehchx1                                    &
         +fuel2(icel)*ehchx2                                    &
         +fuel3(icel)*ehgaze(ico ,i+1)                          &
         +fuel4(icel)*ehgaze(ih2s,i+1)                          &
         +fuel5(icel)*ehgaze(ihy ,i+1)                          &
         +fuel6(icel)*ehgaze(ihcn,i+1)                          &
         +fuel7(icel)*ehgaze(inh3,i+1)                          &
         +oxyd(icel) *ehgaze(io2 ,i+1)                          &
         +prod1(icel)*ehgaze(ico2,i+1)                          &
         +prod2(icel)*ehgaze(ih2o,i+1)                          &
         +prod3(icel)*ehgaze(iso2,i+1)                          &
         +xiner(icel)*ehgaze(in2 ,i+1)

    if ( eh(icel).ge.eh0 .and. eh(icel).le.eh1 ) then
      tp(icel)= th(i) + (eh(icel)-eh0) *                        &
                        (th(i+1)-th(i))/(eh1-eh0)
    endif
  enddo
enddo

!----
! End
!----

return
end subroutine
