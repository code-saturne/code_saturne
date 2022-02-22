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
!> \file cs_fuel_masstransfer.f90
!>
!> \brief Calcultaion of mass transfer terms between the contineous phase
!> and the dispersed phase
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncel          number of cells
!______________________________________________________________________________!

subroutine cs_fuel_masstransfer &
 ( ncel )

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
use cpincl
use ppincl
use ppcpfu
use cs_fuel_incl
use field
!===============================================================================

implicit none

! Arguments

integer          ncel

! Local variables

integer          iel    , icla
integer          ifcvsl

double precision xng,xnuss
double precision pparo2 , xdffli , xdfext , xdftot0 , xdftot1
double precision diacka
double precision dcoke , surf , lambda, visls_0
!
double precision  pref

double precision dhet1, dhet2
double precision deva1, deva2
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cpro_cp, cpro_viscls
double precision, dimension(:), pointer :: cvara_yfolcl, cvara_ngcl
double precision, dimension(:), pointer :: cpro_cgev, cpro_cght, cpro_chgl
double precision, dimension(:), pointer :: cpro_yox, cpro_temp
double precision, dimension(:), pointer :: cpro_diam2, cpro_rom2, cpro_temp2
double precision, dimension(:), pointer :: cpro_rom1

!===============================================================================
! 1. Initializations and preliminary calculations
!===============================================================================
! --- Initialization of mass transfert terms

do icla = 1, nclafu
  call field_get_val_s(igmeva(icla),cpro_cgev)
  call field_get_val_s(igmhtf(icla),cpro_cght)
  call field_get_val_s(ih1hlf(icla),cpro_chgl)
  do iel = 1, ncel
    cpro_cgev(iel) = zero
    cpro_cght(iel) = zero
    cpro_chgl(iel) = zero
  enddo
enddo

! --- Pointer

call field_get_val_s(icrom, crom)
call field_get_val_s(iym1(io2), cpro_yox)
call field_get_val_s(itemp, cpro_temp)
call field_get_val_s(irom1, cpro_rom1)
!
pref = 1.013d0

!===============================================================================
! 2. Source terms for liquid enthalpy
!===============================================================================
!
! Contribution to explicit and implicit balances
! of exchanges by molecular diffusion
! 6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)
!
call field_get_key_int (ivarfl(isca(iscalt)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

do icla = 1, nclafu

  call field_get_val_s(irom2(icla),cpro_rom2)
  call field_get_val_s(idiam2(icla),cpro_diam2)
  call field_get_val_s(itemp2(icla),cpro_temp2)
  call field_get_val_s(igmhtf(icla),cpro_cght)
  call field_get_val_s(ih1hlf(icla),cpro_chgl)

  xnuss = 2.d0

  call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfolcl)
  if ( icp.ge.0 ) call field_get_val_s(icp, cpro_cp)

  call field_get_key_double(ivarfl(isca(iscalt)), kvisl0, visls_0)

  do iel = 1, ncel
    if (ifcvsl.ge.0) then
      if (icp.ge.0) then
        lambda = cpro_viscls(iel) * cpro_cp(iel)
      else
        lambda = cpro_viscls(iel) * cp0
      endif
    else
      if (icp.ge.0) then
        lambda = visls_0 * cpro_cp(iel)
      else
        lambda = visls_0 * cp0
      endif
    endif

    if (       cvara_yfolcl(iel) .gt. epsifl                        &
        .and. cpro_temp(iel).gt. cpro_temp2(iel)) then

       cpro_chgl(iel) = 6.d0*lambda*xnuss/cpro_diam2(iel)**2        &
                           /cpro_rom2(iel)*cvara_yfolcl(iel)

    endif

  enddo

enddo

!===============================================================================
! 3. Fuel mass transfert
!===============================================================================

do icla = 1, nclafu

  call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfolcl)
  call field_get_val_prev_s(ivarfl(isca(ing(icla))), cvara_ngcl)
  call field_get_val_s(itemp2(icla),cpro_temp2)
  call field_get_val_s(igmhtf(icla),cpro_cght)
  call field_get_val_s(ih1hlf(icla),cpro_chgl)
  call field_get_val_s(igmeva(icla),cpro_cgev)

  do iel = 1, ncel

    cpro_cgev(iel) = zero
    cpro_cght(iel) = zero

    if (cvara_yfolcl(iel) .gt. epsifl ) then

!===============================================================================
! Evaporation
!===============================================================================
! Checking of the liquid fuel mass.
! a) If deva1 < deva2 there is no remaining liquid fuel.
!    So there is no evaporation.
! b) Checking of evaporation temperature range.
! c) Checking of the gas phase and drop temperature. It is necessary to have
!     Tgaz > Tgoutte

    deva1 =  cvara_yfolcl(iel)                                             &
             /(cvara_ngcl(iel)*rho0fl)
    deva2 =  (pi*(diniin(icla)**3)/6.d0)+(pi*(dinikf(icla)**3)/6.d0)

      if ( cpro_temp2(iel) .gt. tevap1            .and.                    &
           cpro_temp(iel)  .gt. cpro_temp2(iel) .and. deva1.gt.deva2) then

        ! The evaporated mass flux is determined evapore est determined according
        !               to a supposed profil of dMeva/dTgoutte.

        cpro_cgev(iel) =   cpro_chgl(iel)                                  &
                         / ( hrfvap + cp2fol*(tevap2-cpro_temp2(iel)))

      endif

!===============================================================================
! Heterogeneous combustion
!===============================================================================
! Checking of coke mass.
! a) If deva1.le.deva2 -> There is no remaining liquid fuel. So the particle
!    is composed of coal and inerts.
! b) If dhet1.gt.dhet2 it remains coal. If not the particle is only composed of
!    inerts.

    dhet1= cvara_yfolcl(iel)                                               &
            /(cvara_ngcl(iel)*rho0fl)
    dhet2= pi*(diniin(icla)**3)/6.d0

      if (deva1.le.deva2.and.dhet1.gt.dhet2 ) then

      ! We consider that the remaining coke mass as a spheric particle. The
      !   correspondent diameter is dcoke.
      dcoke = ( ( cvara_yfolcl(iel)                                        &
              /(cvara_ngcl(iel)*rho0fl)                                    &
              -pi*(diniin(icla)**3)/6.d0  )                                &
               *6.d0/pi )**(1.d0/3.d0)

      ! Calculation of the partial pressure of oxygene [atm]                                                 ---
      !   PO2 = RHO1*CS_PHYSICAL_CONSTANTS_R*T*YO2/MO2

      pparo2 = cpro_rom1(iel)*cs_physical_constants_r                      &
              *cpro_temp(iel)*cpro_yox(iel)/wmole(io2)
      pparo2 = pparo2 / prefth

      ! Chemical kinetic coefficient of CO forming
      !   in [kg.m-2.s-1.atm(-n)]
      xdffli = ahetfl*exp(-ehetfl*4185.d0                                  &
              /(cs_physical_constants_r*cpro_temp(iel)))

      ! Coefficient of diffusion in kg/m2/s/[atm]: XDFEXT
      ! Global coefficient for n=0.5 in kg/m2/s: XDFTOT0
      ! Global coefficient for n=1 in Kg/m2/s: XDFTOT1

      diacka = dcoke/(dinikf(icla))
        if ( diacka .gt. epsifl ) then
        xdfext = 2.53d-7*((cpro_temp(iel))**0.75d0) / dcoke*2.d0
        xdftot1 = pparo2 / ( 1.d0/xdffli + 1.d0/xdfext )
        xdftot0 = -(xdffli**2)/(2.d0*xdfext**2)+(pparo2*xdffli**2          &
                +(xdffli**4)/(2.d0*xdfext**2))**0.5d0
        else
        xdftot1 = xdffli*pparo2
        xdftot0 = xdffli*pparo2**0.5d0
        endif

      ! Spherical particle surface.
      surf = pi*(dcoke**2)

      ! Number of particles in the cell.
      xng   = cvara_ngcl(iel)

        if (iofhet.eq.1) then
          cpro_cght(iel) = - xdftot1*surf*xng
        else
          cpro_cght(iel) = - xdftot0*surf*xng
        endif

      endif

    endif

  enddo

enddo

!----
! End
!----

return
end subroutine
