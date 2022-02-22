!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!> \file cs_coal_thfieldconv1.f90
!>
!> \brief Calculation of the gas temperature
!>        Function with gas enthalpy and concentrations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     location_id     mesh location id (cells or boundary faces)
!> \param[in]     eh              gas enthalpy
!>                                (j/kg of gaseous mixture)
!> \param[in,out] tp              gas temperature in kelvin
!______________________________________________________________________________!

subroutine cs_coal_thfieldconv1(location_id, eh, tp)  &
 bind(C, name='cs_coal_thfieldconv1')

!==============================================================================
! Module files
!==============================================================================

use, intrinsic :: iso_c_binding

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
use mesh
use field
use cs_c_bindings
use pointe

!===============================================================================

implicit none

! Arguments

integer(c_int), value :: location_id
real(kind=c_double), dimension(*) :: eh, tp

! Local variables

integer          i, iel, ielt, nelt, icha

double precision ychx10 , ychx20 , ehchx1 , ehchx2
double precision den1   , den2 , eh0 , eh1
double precision f1mc(ncharm), f2mc(ncharm)

double precision, dimension(:), pointer :: x1
double precision, dimension(:), pointer :: fuel1, fuel2, fuel3, fuel4, fuel5
double precision, dimension(:), pointer :: fuel6, fuel7, oxyd
double precision, dimension(:), pointer :: prod1, prod2, prod3 , xiner
type(pmapper_double_r1), dimension(:), allocatable :: cvar_f1m, cvar_f2m

!===============================================================================

call field_get_val_s(iym1(ichx1), fuel1)
call field_get_val_s(iym1(ichx2), fuel2)
call field_get_val_s(iym1(ico  ), fuel3)
call field_get_val_s(iym1(ih2s ), fuel4)
call field_get_val_s(iym1(ihy  ), fuel5)
call field_get_val_s(iym1(ihcn ), fuel6)
call field_get_val_s(iym1(inh3 ), fuel7)
call field_get_val_s(iym1(io2  ), oxyd)
call field_get_val_s(iym1(ico2 ), prod1)
call field_get_val_s(iym1(ih2o ), prod2)
call field_get_val_s(iym1(iso2 ), prod3)
call field_get_val_s(iym1(in2  ), xiner)

! Massic fraction of gas
call field_get_val_s_by_name("x_c", x1)

if (location_id .eq. MESH_LOCATION_CELLS) then
  nelt = ncel
else if (location_id .eq. MESH_LOCATION_BOUNDARY_FACES) then
  nelt = nfabor
else
  call csexit(1)
endif

!===============================================================================

allocate(cvar_f1m(ncharb))
allocate(cvar_f2m(ncharb))
do icha = 1, ncharb
  call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m(icha)%p)
  call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m(icha)%p)
enddo

i = npo-1

do ielt = 1, nelt

  if (location_id .eq. MESH_LOCATION_CELLS) then
    iel = ielt
  else
    iel = ifabor(ielt)
  endif

! --- Calculation of enthalpy of the gaseous species CHx1m
!                                            and CHx2m at TH(NPO)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb

    f1mc(icha) = cvar_f1m(icha)%p(iel) / x1(iel)
    f2mc(icha) = cvar_f2m(icha)%p(iel) / x1(iel)

    den1   = 1.d0                                                  &
         / ( a1(icha)*wmole(ichx1c(icha))                          &
            +b1(icha)*wmole(ico)                                   &
            +c1(icha)*wmole(ih2o)                                  &
            +d1(icha)*wmole(ih2s)                                  &
            +e1(icha)*wmole(ihcn)                                  &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10                                                &
            +den1*(f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1                                                &
            +den1*( ehgaze(ichx1c(icha),i+1)                       &
                   *f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                  &
         / ( a2(icha)*wmole(ichx2c(icha))                          &
            +b2(icha)*wmole(ico)                                   &
            +c2(icha)*wmole(ih2o)                                  &
            +d2(icha)*wmole(ih2s)                                  &
            +e2(icha)*wmole(ihcn)                                  &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20                                                &
            +den2 *(f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                       &
         ( ehgaze(ichx2c(icha),i+1)                                &
          *f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if (ychx10.gt.epzero) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i+1)
  endif
  if (ychx20.gt.epzero) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i+1)
  endif

  ! --- Eventual clipping of temperature at TH(NPO) if EH > EH1

  eh1 = fuel1(iel)*ehchx1                                  &
       +fuel2(iel)*ehchx2                                  &
       +fuel3(iel)*ehgaze(ico ,i+1)                        &
       +fuel4(iel)*ehgaze(ih2s,i+1)                        &
       +fuel5(iel)*ehgaze(ihy ,i+1)                        &
       +fuel6(iel)*ehgaze(ihcn,i+1)                        &
       +fuel7(iel)*ehgaze(inh3,i+1)                        &
       +oxyd(iel) *ehgaze(io2 ,i+1)                        &
       +prod1(iel)*ehgaze(ico2,i+1)                        &
       +prod2(iel)*ehgaze(ih2o,i+1)                        &
       +prod3(iel)*ehgaze(iso2,i+1)                        &
       +xiner(iel)*ehgaze(in2 ,i+1)

  if (eh(ielt).ge.eh1) tp(ielt) = th(i+1)
enddo

i = 1

do ielt = 1, nelt

  if (location_id .eq. MESH_LOCATION_CELLS) then
    iel = ielt
  else
    iel = ifabor(ielt)
  endif

  ! --- Calculation of enthalpy of the gaseous species CHx1m
  !                                            and CHx2m at TH(1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    f1mc(icha) = cvar_f1m(icha)%p(iel) / x1(iel)
    f2mc(icha) = cvar_f2m(icha)%p(iel) / x1(iel)
    den1   = 1.d0                                                    &
         / ( a1(icha)*wmole(ichx1c(icha))                            &
            +b1(icha)*wmole(ico)                                     &
            +c1(icha)*wmole(ih2o)                                    &
            +d1(icha)*wmole(ih2s)                                    &
            +e1(icha)*wmole(ihcn)                                    &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10                                                  &
            +den1*(f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1                                                  &
            +den1*( ehgaze(ichx1c(icha),i)                           &
                   *f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                     &
         / ( a2(icha)*wmole(ichx2c(icha))                             &
            +b2(icha)*wmole(ico)                                      &
            +c2(icha)*wmole(ih2o)                                     &
            +d2(icha)*wmole(ih2s)                                     &
            +e2(icha)*wmole(ihcn)                                     &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20                                                   &
            +den2*(f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2                                                   &
            +den2*( ehgaze(ichx2c(icha),i)                            &
                   *f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
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

  eh0  = fuel1(iel)*ehchx1                                  &
        +fuel2(iel)*ehchx2                                  &
        +fuel3(iel)*ehgaze(ico ,i)                          &
        +fuel4(iel)*ehgaze(ih2s,i)                          &
        +fuel5(iel)*ehgaze(ihy ,i)                          &
        +fuel6(iel)*ehgaze(ihcn,i)                          &
        +fuel7(iel)*ehgaze(inh3,i)                          &
        +oxyd(iel) *ehgaze(io2 ,i)                          &
        +prod1(iel)*ehgaze(ico2,i)                          &
        +prod2(iel)*ehgaze(ih2o,i)                          &
        +prod3(iel)*ehgaze(iso2,i)                          &
        +xiner(iel)*ehgaze(in2 ,i)

  if (eh(ielt) .le. eh0 ) tp(ielt)= th(i)

enddo


do i = 1, npo-1

  do ielt = 1, nelt

    if (location_id .eq. MESH_LOCATION_CELLS) then
      iel = ielt
    else
      iel = ifabor(ielt)
    endif

    ! --- Calculation of enthalpy of the gaseous species CHx1m
    !                                            and CHx2m for TH(I)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      f1mc(icha) = cvar_f1m(icha)%p(iel) / x1(iel)
      f2mc(icha) = cvar_f2m(icha)%p(iel) / x1(iel)
    enddo
    do icha = 1, ncharb
      den1   = 1.d0                                                 &
             / ( a1(icha)*wmole(ichx1c(icha))                       &
                +b1(icha)*wmole(ico)                                &
                +c1(icha)*wmole(ih2o)                               &
                +d1(icha)*wmole(ih2s)                               &
                +e1(icha)*wmole(ihcn)                               &
                +f1(icha)*wmole(inh3) )
      ychx10 = ychx10                                               &
              +den1*f1mc(icha)*a1(icha)*wmole(ichx1c(icha))
      ehchx1 = ehchx1                                               &
              +den1*( ehgaze(ichx1c(icha),i)                        &
                     *f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                                 &
             / ( a2(icha)*wmole(ichx2c(icha))                       &
                +b2(icha)*wmole(ico)                                &
                +c2(icha)*wmole(ih2o)                               &
                +d2(icha)*wmole(ih2s)                               &
                +e2(icha)*wmole(ihcn)                               &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20                                               &
              +den2*(f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2                                               &
              +den2*( ehgaze(ichx2c(icha),i)                        &
                     *f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
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
    eh0 = fuel1(iel)*ehchx1                                  &
         +fuel2(iel)*ehchx2                                  &
         +fuel3(iel)*ehgaze(ico ,i)                          &
         +fuel4(iel)*ehgaze(ih2s,i)                          &
         +fuel5(iel)*ehgaze(ihy ,i)                          &
         +fuel6(iel)*ehgaze(ihcn,i)                          &
         +fuel7(iel)*ehgaze(inh3,i)                          &
         +oxyd(iel) *ehgaze(io2 ,i)                          &
         +prod1(iel)*ehgaze(ico2,i)                          &
         +prod2(iel)*ehgaze(ih2o,i)                          &
         +prod3(iel)*ehgaze(iso2,i)                          &
         +xiner(iel)*ehgaze(in2 ,i)

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
              +den1*f1mc(icha)*a1(icha)*wmole(ichx1c(icha))
      ehchx1 = ehchx1                                                 &
              +den1*( ehgaze(ichx1c(icha),i+1)                        &
                     *f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                                   &
             / ( a2(icha)*wmole(ichx2c(icha))                         &
                +b2(icha)*wmole(ico)                                  &
                +c2(icha)*wmole(ih2o)                                 &
                +d2(icha)*wmole(ih2s)                                 &
                +e2(icha)*wmole(ihcn)                                 &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20                                                 &
              +den2*f2mc(icha)*a2(icha)*wmole(ichx2c(icha))
      ehchx2 = ehchx2                                                 &
              +den2*( ehgaze(ichx2c(icha),i+1)                        &
                     *f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if (ychx10.gt.epzero) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i+1)
    endif
    if (ychx20.gt.epzero) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i+1)
    endif

    eh1 = fuel1(iel)*ehchx1                                    &
         +fuel2(iel)*ehchx2                                    &
         +fuel3(iel)*ehgaze(ico ,i+1)                          &
         +fuel4(iel)*ehgaze(ih2s,i+1)                          &
         +fuel5(iel)*ehgaze(ihy ,i+1)                          &
         +fuel6(iel)*ehgaze(ihcn,i+1)                          &
         +fuel7(iel)*ehgaze(inh3,i+1)                          &
         +oxyd(iel) *ehgaze(io2 ,i+1)                          &
         +prod1(iel)*ehgaze(ico2,i+1)                          &
         +prod2(iel)*ehgaze(ih2o,i+1)                          &
         +prod3(iel)*ehgaze(iso2,i+1)                          &
         +xiner(iel)*ehgaze(in2 ,i+1)

    if (eh(ielt).ge.eh0 .and. eh(ielt).le.eh1) then
      tp(ielt)= th(i) + (  eh(ielt)-eh0)                       &
                         * (th(i+1)-th(i))/(eh1-eh0)
    endif
  enddo
enddo

!----
! End
!----

return
end subroutine
