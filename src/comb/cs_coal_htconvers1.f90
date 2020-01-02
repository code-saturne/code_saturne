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
! Function :
! --------
!> \file cs_coal_htconvers1.f90
!>
!> \brief - Calculation of gas temperature
!>          Function with gas enthalpy and concentrations
!>          if mode = 1
!>        - Calculation of gas enthalpy
!>          Function with gas temperature and concentrations
!>          if mode = -1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     mode          -1 : t -> h  ;   1 : h -> t
!> \param[in,out] eh            gas enthalpy
!>                                   (\f$ j . kg^{-1} \f$ of mixed gas)
!> \param[in]     xesp          mass fraction of species
!> \param[in]     f1mc          average f1
!> \param[in]     f2mc          average f2
!> \param[in,out] tp            gas temperature (in kelvin)
!______________________________________________________________________________!

subroutine cs_coal_htconvers1 &
 ( mode , eh , xesp , f1mc , f2mc , tp )

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
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          mode
double precision eh,tp
double precision xesp(ngazem)
double precision f1mc(ncharm),f2mc(ncharm)
! Local variables

integer          i , icha

double precision ychx10 , ychx20 , ehchx1 , ehchx2
double precision den1   , den2
double precision eh0 , eh1

!===============================================================================
! 1. Calculation of temperature from enthalpy
!===============================================================================

if (mode .eq. 1) then

  i = npo

  ! Calculation of enthalpy of the gaseous species CHx1m
  !                                            and CHx2m at TH(NPO)

  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero

  do icha = 1, ncharb
    den1   = 1.d0                                                 &
         / ( a1(icha)*wmole(ichx1c(icha))                         &
            +b1(icha)*wmole(ico )                                 &
            +c1(icha)*wmole(ih2o)                                 &
            +d1(icha)*wmole(ih2s)                                 &
            +e1(icha)*wmole(ihcn)                                 &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10 + den1 *                                      &
         ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
         ( ehgaze(ichx1c(icha),i)*                                &
         f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
         / ( a2(icha)*wmole(ichx2c(icha))                         &
            +b2(icha)*wmole(ico )                                 &
            +c2(icha)*wmole(ih2o)                                 &
            +d2(icha)*wmole(ih2s)                                 &
            +e2(icha)*wmole(ihcn)                                 &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20 + den2 *                                      &
         ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
         ( ehgaze(ichx2c(icha),i)*                                &
         f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo

  if (ychx10.gt.epzero) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i)
  endif
  if (ychx20.gt.epzero) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i)
  endif

! --- Eventual clipping of temperature at TH(NPO) if EH > EH1
  eh1 = xesp(ichx1)*ehchx1                                         &
       + xesp(ichx2)*ehchx2                                        &
       + xesp(ico  )*ehgaze(ico  ,i)                               &
       + xesp(ih2s )*ehgaze(ih2s ,i)                               &
       + xesp(ihy  )*ehgaze(ihy  ,i)                               &
       + xesp(ihcn )*ehgaze(ihcn ,i)                               &
       + xesp(io2  )*ehgaze(io2  ,i)                               &
       + xesp(ico2 )*ehgaze(ico2 ,i)                               &
       + xesp(ih2o )*ehgaze(ih2o ,i)                               &
       + xesp(iso2 )*ehgaze(iso2 ,i)                               &
       + xesp(inh3 )*ehgaze(inh3 ,i)                               &
       + xesp(in2  )*ehgaze(in2  ,i)

  if (eh .ge. eh1) then
    tp = th(i)
    goto 501
  endif

  i = 1

  ! Calculation of enthalpy of the gaseous species CHx1m
  !                                             and CHx2m at TH(1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero

  do icha = 1, ncharb
    den1   = 1.d0                                                      &
         / ( a1(icha)*wmole(ichx1c(icha))                              &
            +b1(icha)*wmole(ico )                                      &
            +c1(icha)*wmole(ih2o)                                      &
            +d1(icha)*wmole(ih2s)                                      &
            +e1(icha)*wmole(ihcn)                                      &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10 + den1*(f1mc(icha)*a1(icha)*wmole(ichx1c(icha)))
    ehchx1 = ehchx1 + den1*(   ehgaze(ichx1c(icha),i)                  &
                            * f1mc(icha)*a1(icha)*wmole(ichx1c(icha)))
    den2   = 1.d0                                                      &
         / ( a2(icha)*wmole(ichx2c(icha))                         &
            +b2(icha)*wmole(ico )                                 &
            +c2(icha)*wmole(ih2o)                                 &
            +d2(icha)*wmole(ih2s)                                 &
            +e2(icha)*wmole(ihcn)                                 &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20 + den2*(f2mc(icha)*a2(icha)*wmole(ichx2c(icha)))
    ehchx2 = ehchx2 + den2*(  ehgaze(ichx2c(icha),i)                   &
                            * f2mc(icha)*a2(icha)*wmole(ichx2c(icha)))
  enddo
  if (ychx10.gt.epzero) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i)
  endif
  if (ychx20.gt.epzero) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i)
  endif

  ! --- Eventual clipping of temperature at TH(1) if EH < EH0
  eh0 = xesp(ichx1)*ehchx1                                        &
      + xesp(ichx2)*ehchx2                                        &
      + xesp(ico  )*ehgaze(ico  ,i)                               &
      + xesp(ih2s )*ehgaze(ih2s ,i)                               &
      + xesp(ihy  )*ehgaze(ihy  ,i)                               &
      + xesp(ihcn )*ehgaze(ihcn ,i)                               &
      + xesp(io2  )*ehgaze(io2  ,i)                               &
      + xesp(ico2 )*ehgaze(ico2 ,i)                               &
      + xesp(ih2o )*ehgaze(ih2o ,i)                               &
      + xesp(iso2 )*ehgaze(iso2 ,i)                               &
      + xesp(inh3 )*ehgaze(inh3 ,i)                               &
      + xesp(in2  )*ehgaze(in2  ,i)

  if ( eh .le. eh0 ) then
    tp= th(i)
    goto 501
  endif

 500    continue
  i = i + 1

! --- Calculation of enthalpy of the gaseous species CHx1m
!                                            and CHx2m for TH(I-1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero

  do icha = 1, ncharb
    den1   = 1.d0                                                       &
         / ( a1(icha)*wmole(ichx1c(icha))                               &
            +b1(icha)*wmole(ico )                                       &
            +c1(icha)*wmole(ih2o)                                       &
            +d1(icha)*wmole(ih2s)                                       &
            +e1(icha)*wmole(ihcn)                                       &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10 + den1*(f1mc(icha)*a1(icha)*wmole(ichx1c(icha)))
    ehchx1 = ehchx1 + den1*(   ehgaze(ichx1c(icha),i-1)                 &
                            * f1mc(icha)*a1(icha)*wmole(ichx1c(icha)))
    den2   = 1.d0                                                       &
          / ( a2(icha)*wmole(ichx2c(icha))                              &
             +b2(icha)*wmole(ico )                                      &
             +c2(icha)*wmole(ih2o)                                      &
             +d2(icha)*wmole(ih2s)                                      &
             +e2(icha)*wmole(ihcn)                                      &
             +f2(icha)*wmole(inh3) )
    ychx20 = ychx20 + den2*(f2mc(icha)*a2(icha)*wmole(ichx2c(icha)))
    ehchx2 = ehchx2 + den2*(  ehgaze(ichx2c(icha),i-1)                  &
                            * f2mc(icha)*a2(icha)*wmole(ichx2c(icha)))
  enddo

  if (ychx10.gt.epzero) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i-1)
  endif
  if (ychx20.gt.epzero) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i-1)
  endif

  eh0 = xesp(ichx1)*ehchx1                                          &
      + xesp(ichx2)*ehchx2                                          &
      + xesp(ico  )*ehgaze(ico  ,i-1)                               &
      + xesp(ih2s )*ehgaze(ih2s ,i-1)                               &
      + xesp(ihy  )*ehgaze(ihy  ,i-1)                               &
      + xesp(ihcn )*ehgaze(ihcn ,i-1)                               &
      + xesp(io2  )*ehgaze(io2  ,i-1)                               &
      + xesp(ico2 )*ehgaze(ico2 ,i-1)                               &
      + xesp(ih2o )*ehgaze(ih2o ,i-1)                               &
      + xesp(iso2 )*ehgaze(iso2 ,i-1)                               &
      + xesp(inh3 )*ehgaze(inh3 ,i-1)                               &
      + xesp(in2  )*ehgaze(in2  ,i-1)

  ! Calculation of the enthalpy of the gaseous species CHx1m
  !                                                and CHx2m at TH(I)

  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
         / ( a1(icha)*wmole(ichx1c(icha))                         &
            +b1(icha)*wmole(ico )                                 &
            +c1(icha)*wmole(ih2o)                                 &
            +d1(icha)*wmole(ih2s)                                 &
            +e1(icha)*wmole(ihcn)                                 &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10 + den1 *                                      &
      ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
      ( ehgaze(ichx1c(icha),i)*                                   &
        f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
         / ( a2(icha)*wmole(ichx2c(icha))                         &
            +b2(icha)*wmole(ico )                                 &
            +c2(icha)*wmole(ih2o)                                 &
            +d2(icha)*wmole(ih2s)                                 &
            +e2(icha)*wmole(ihcn)                                 &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20 + den2 *                                      &
       ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
       ( ehgaze(ichx2c(icha),i)*                                  &
         f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
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

  eh1 = xesp(ichx1)*ehchx1                                        &
      + xesp(ichx2)*ehchx2                                        &
      + xesp(ico  )*ehgaze(ico  ,i)                               &
      + xesp(ih2s )*ehgaze(ih2s ,i)                               &
      + xesp(ihy  )*ehgaze(ihy  ,i)                               &
      + xesp(ihcn )*ehgaze(ihcn ,i)                               &
      + xesp(io2  )*ehgaze(io2  ,i)                               &
      + xesp(ico2 )*ehgaze(ico2 ,i)                               &
      + xesp(ih2o )*ehgaze(ih2o ,i)                               &
      + xesp(iso2 )*ehgaze(iso2 ,i)                               &
      + xesp(inh3 )*ehgaze(inh3 ,i)                               &
      + xesp(in2  )*ehgaze(in2  ,i)

  if ( eh .ge. eh0  .and. eh .le. eh1  ) then
    tp = th(i-1) + (eh-eh0) *                                     &
                 (th(i)-th(i-1))/(eh1-eh0)
    goto 501
  endif
  goto 500

 501    continue

!===============================================================================
! 2. Calculation of the enthalpy from temperature
!===============================================================================

else if ( mode .eq. -1 ) then

  i = npo
  ! --- Calculation at Max
  if ( tp .ge. th(i) ) then
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
           / ( a1(icha)*wmole(ichx1c(icha))                       &
              +b1(icha)*wmole(ico)                                &
              +c1(icha)*wmole(ih2o)                               &
              +d1(icha)*wmole(ih2s)                               &
              +e1(icha)*wmole(ihcn)                               &
              +f1(icha)*wmole(inh3) )
      ychx10 = ychx10 + den1 *                                    &
           ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
           ( ehgaze(ichx1c(icha),i)*                              &
           f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
           / ( a2(icha)*wmole(ichx2c(icha))                       &
              +b2(icha)*wmole(ico)                                &
              +c2(icha)*wmole(ih2o)                               &
              +d2(icha)*wmole(ih2s)                               &
              +e2(icha)*wmole(ihcn)                               &
              +f2(icha)*wmole(inh3) )
      ychx20 = ychx20 + den2 *                                    &
           ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
           ( ehgaze(ichx2c(icha),i)*                              &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
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
    eh = xesp(ichx1)*ehchx1                                       &
       + xesp(ichx2)*ehchx2                                       &
       + xesp(ico  )*ehgaze(ico  ,i)                               &
       + xesp(ih2s )*ehgaze(ih2s ,i)                               &
       + xesp(ihy  )*ehgaze(ihy  ,i)                               &
       + xesp(ihcn )*ehgaze(ihcn ,i)                               &
       + xesp(io2  )*ehgaze(io2  ,i)                               &
       + xesp(ico2 )*ehgaze(ico2 ,i)                               &
       + xesp(ih2o )*ehgaze(ih2o ,i)                               &
       + xesp(iso2 )*ehgaze(iso2 ,i)                               &
       + xesp(inh3 )*ehgaze(inh3 ,i)                               &
       + xesp(in2  )*ehgaze(in2  ,i)
    goto 601
  endif
  ! Clipping at Min
  i = 1
  if ( tp .le. th(i) ) then
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
           / ( a1(icha)*wmole(ichx1c(icha))                       &
              +b1(icha)*wmole(ico )                               &
              +c1(icha)*wmole(ih2o)                               &
              +d1(icha)*wmole(ih2s)                               &
              +e1(icha)*wmole(ihcn)                               &
              +f1(icha)*wmole(inh3) )
      ychx10 = ychx10 + den1 *                                    &
           ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
           ( ehgaze(ichx1c(icha),i)*                              &
           f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
           / ( a2(icha)*wmole(ichx2c(icha))                       &
              +b2(icha)*wmole(ico )                               &
              +c2(icha)*wmole(ih2o)                               &
              +d2(icha)*wmole(ih2s)                               &
              +e2(icha)*wmole(ihcn)                               &
              +f2(icha)*wmole(inh3 ) )
      ychx20 = ychx20 + den2 *                                    &
           ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
           ( ehgaze(ichx2c(icha),i)*                              &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
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
    eh = xesp(ichx1)*ehchx1                                       &
       + xesp(ichx2)*ehchx2                                       &
       + xesp(ico  )*ehgaze(ico  ,i)                               &
       + xesp(ih2s )*ehgaze(ih2s ,i)                               &
       + xesp(ihy  )*ehgaze(ihy  ,i)                               &
       + xesp(ihcn )*ehgaze(ihcn ,i)                               &
       + xesp(io2  )*ehgaze(io2  ,i)                               &
       + xesp(ico2 )*ehgaze(ico2 ,i)                               &
       + xesp(ih2o )*ehgaze(ih2o ,i)                               &
       + xesp(iso2 )*ehgaze(iso2 ,i)                               &
       + xesp(inh3 )*ehgaze(inh3 ,i)                               &
       + xesp(in2  )*ehgaze(in2  ,i)
    goto 601
  endif
 ! Interpolation in the table
  i = 1
 600    continue

  i = i + 1
  if ( tp .le. th(i) ) then
    ! --- Calculation of the enthalpy of the gaseous species TH(I-1)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
             / ( a1(icha)*wmole(ichx1c(icha))                     &
              +b1(icha)*wmole(ico)                                &
              +c1(icha)*wmole(ih2o)                               &
              +d1(icha)*wmole(ih2s)                               &
              +e1(icha)*wmole(ihcn)                               &
              +f1(icha)*wmole(inh3) )
      ychx10 = ychx10 + den1 *                                    &
        ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
        ( ehgaze(ichx1c(icha),i-1)*                               &
          f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
             / ( a2(icha)*wmole(ichx2c(icha))                     &
                +b2(icha)*wmole(ico)                              &
                +c2(icha)*wmole(ih2o)                             &
                +d2(icha)*wmole(ih2s)                             &
                +e2(icha)*wmole(ihcn)                             &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20 + den2 *                                    &
         ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
         ( ehgaze(ichx2c(icha),i-1)*                              &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i-1)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i-1)
    endif

    eh0 = xesp(ichx1)*ehchx1                                      &
        + xesp(ichx2)*ehchx2                                      &
        + xesp(ico  )*ehgaze(ico  ,i-1)                           &
        + xesp(ih2s )*ehgaze(ih2s ,i-1)                           &
        + xesp(ihy  )*ehgaze(ihy  ,i-1)                           &
        + xesp(ihcn )*ehgaze(ihcn ,i-1)                           &
        + xesp(io2  )*ehgaze(io2  ,i-1)                           &
        + xesp(ico2 )*ehgaze(ico2 ,i-1)                           &
        + xesp(ih2o )*ehgaze(ih2o ,i-1)                           &
        + xesp(iso2 )*ehgaze(iso2 ,i-1)                           &
        + xesp(inh3 )*ehgaze(inh3 ,i-1)                           &
        + xesp(in2  )*ehgaze(in2  ,i-1)
    ! --- Calculation of the enthalpy of the gaseous species TH(I)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                              &
             / ( a1(icha)*wmole(ichx1c(icha))                    &
                +b1(icha)*wmole(ico)                             &
                +c1(icha)*wmole(ih2o)                            &
                +d1(icha)*wmole(ih2s)                            &
                +e1(icha)*wmole(ihcn)                            &
                +f1(icha)*wmole(inh3) )
      ychx10 = ychx10 + den1 *                                   &
        ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                   &
        ( ehgaze(ichx1c(icha),i)*                                &
          f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                              &
             / ( a2(icha)*wmole(ichx2c(icha))                    &
                +b2(icha)*wmole(ico)                             &
                +c2(icha)*wmole(ih2o)                            &
                +d2(icha)*wmole(ih2s)                            &
                +e2(icha)*wmole(ihcn)                            &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20 + den2 *                                   &
         ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                   &
         ( ehgaze(ichx2c(icha),i)*                               &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
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

    eh1 = xesp(ichx1)*ehchx1                                     &
        + xesp(ichx2)*ehchx2                                     &
        + xesp(ico  )*ehgaze(ico  ,i)                            &
        + xesp(ih2s )*ehgaze(ih2s ,i)                            &
        + xesp(ihy  )*ehgaze(ihy  ,i)                            &
        + xesp(ihcn )*ehgaze(ihcn ,i)                            &
        + xesp(io2  )*ehgaze(io2  ,i)                            &
        + xesp(ico2 )*ehgaze(ico2 ,i)                            &
        + xesp(ih2o )*ehgaze(ih2o ,i)                            &
        + xesp(iso2 )*ehgaze(iso2 ,i)                            &
        + xesp(inh3 )*ehgaze(inh3 ,i)                            &
        + xesp(in2  )*ehgaze(in2  ,i)

    eh =  eh0                                                    &
           + (eh1-eh0)*(tp-th(i-1))/(th(i)-th(i-1))
    goto 601
  endif
  goto 600

 601    continue

else
  write(nfecra,1000) mode
  call csexit(1)
endif

!--------
! Formats
!--------
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Error in cpthp1                                ',/,&
'@    =========                                               ',/,&
'@    Incorrect value of argument mode                        ',/,&
'@    it must be an integer equal to 1 or -1                  ',/,&
'@    it worths here ',I10                                     ,/,&
'@                                                            ',/,&
'@  The calculation can not run.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
!----
! End
!----

return
end subroutine
