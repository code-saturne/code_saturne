!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in,out] propce        physic properties at cell centers
!______________________________________________________________________________!

subroutine cs_fuel_masstransfer &
 ( ncelet , ncel   ,            &
   propce )

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

integer          ncelet , ncel

double precision propce(ncelet,*)

! Local variables

integer          iel    , icla
integer          ipcte1 , ipcte2 , ipcro2 , ipcdia
integer          ipcgev , ipcght , ipcyox
integer          ifcvsl , ipchgl

double precision xng,xnuss
double precision pparo2 , xdffli , xdfext , xdftot0 , xdftot1
double precision diacka
double precision dcoke , surf , lambda
!
double precision  pref

double precision dhet1, dhet2
double precision deva1, deva2
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cpro_cp, cpro_viscls
double precision, dimension(:), pointer :: cvara_yfolcl, cvara_ngcl

!===============================================================================
! 1. Initializations and preliminary calculations
!===============================================================================
! --- Initialization of mass transfert terms

do icla = 1, nclafu
  ipcgev = ipproc(igmeva(icla))
  ipcght = ipproc(igmhtf(icla))
  ipchgl = ipproc(ih1hlf(icla))
  do iel = 1, ncel
    propce(iel,ipcgev) = zero
    propce(iel,ipcght) = zero
    propce(iel,ipchgl) = zero
  enddo
enddo

! --- Pointer

call field_get_val_s(icrom, crom)
ipcte1 = ipproc(itemp1)
ipcyox = ipproc(iym1(io2))
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

if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

do icla = 1, nclafu

  ipcro2 = ipproc(irom2 (icla))
  ipcdia = ipproc(idiam2(icla))
  ipcte2 = ipproc(itemp2(icla))
  ipcght = ipproc(igmhtf(icla))
  ipchgl = ipproc(ih1hlf(icla))

  xnuss = 2.d0

  call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfolcl)
  if ( icp.gt.0 ) call field_get_val_s(iprpfl(icp), cpro_cp)

  do iel = 1, ncel
    if (ifcvsl.ge.0) then
      if (icp.gt.0) then
        lambda = cpro_viscls(iel) * cpro_cp(iel)
      else
        lambda = cpro_viscls(iel) * cp0
      endif
    else
      if (icp.gt.0) then
        lambda = visls0(iscalt) * cpro_cp(iel)
      else
        lambda = visls0(iscalt) * cp0
      endif
    endif

    if ( cvara_yfolcl(iel) .gt. epsifl  .and.                              &
         propce(iel,ipcte1).gt. propce(iel,ipcte2)        ) then

       propce(iel,ipchgl) = 6.d0*lambda*xnuss/propce(iel,ipcdia)**2        &
                           /propce(iel,ipcro2)*cvara_yfolcl(iel)

    endif

  enddo

enddo

!===============================================================================
! 3. Fuel mass transfert
!===============================================================================

do icla = 1, nclafu

  call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfolcl)
  call field_get_val_prev_s(ivarfl(isca(ing(icla))), cvara_ngcl)
  ipcro2 = ipproc(irom2 (icla))
  ipcdia = ipproc(idiam2(icla))
  ipcte2 = ipproc(itemp2(icla))
  ipcgev = ipproc(igmeva(icla))
  ipchgl = ipproc(ih1hlf(icla))
  ipcght = ipproc(igmhtf(icla))

  do iel = 1, ncel

    propce(iel,ipcgev) = zero
    propce(iel,ipcght) = zero

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

      if ( propce(iel,ipcte2)    .gt. tevap1               .and.           &
           propce(iel,ipcte1)    .gt. propce(iel,ipcte2)   .and.           &
           deva1.gt.deva2                                          ) then

      ! The evaporated mass flux is determined evapore est determined according
      !               to a supposed profil of dMeva/dTgoutte.

      propce(iel,ipcgev) = propce(iel,ipchgl)                              &
                          /( hrfvap + cp2fol*(tevap2-propce(iel,ipcte2)) )

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
      !   PO2 = RHO1*RR*T*YO2/MO2

      pparo2 = propce(iel,ipproc(irom1))*rr*propce(iel,ipcte1)             &
              *propce(iel,ipcyox)/wmole(io2)
      pparo2 = pparo2 / prefth

      ! Chemical kinetic coefficient of CO forming
      !   in [kg.m-2.s-1.atm(-n)]
      xdffli = ahetfl*exp(-ehetfl*4185.d0                                  &
              /(rr*propce(iel,ipcte1)))

      ! Coefficient of diffusion in kg/m2/s/[atm]: XDFEXT
      ! Global coefficient for n=0.5 in kg/m2/s: XDFTOT0
      ! Global coefficient for n=1 in Kg/m2/s: XDFTOT1

      diacka = dcoke/(dinikf(icla))
        if ( diacka .gt. epsifl ) then
        xdfext = 2.53d-7*((propce(iel,ipcte1))**0.75d0)                    &
                / dcoke*2.d0
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
        propce(iel,ipcght) = - xdftot1*surf*xng
        else
        propce(iel,ipcght) = - xdftot0*surf*xng
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
