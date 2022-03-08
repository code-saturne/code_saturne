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

!===============================================================================
!  Purpose:
!  --------

!> \file kinrates.f90
!> \brief Calls the computation of reaction rates for atmospheric chemistry

!-------------------------------------------------------------------------------
subroutine kinrates ()

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use atincl
use atchem
use cs_c_bindings

implicit none

!===============================================================================

! Local variables

integer iel,ii,iphotolysis
integer          met_qv_id
double precision temp, dens          ! temperature, density
double precision press, hspec        ! pressure, specific humidity (kg/kg)
double precision rk(nrg)             ! kinetic rates
double precision azi                 ! zenith angle
double precision dlmuzero            ! cos of zenith angle
double precision zent                ! Z coordinate of a cell
double precision qureel              ! julian day
double precision heurtu              ! yime (UTC)
double precision albe                ! albedo, useless here
double precision fo                  ! solar constant, useless here
double precision omega

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_totwt
double precision, dimension(:), pointer :: cpro_tempc, cpro_liqwt
double precision, dimension(:), pointer :: cpro_met_temp
double precision, dimension(:), pointer :: cpro_met_p
double precision, dimension(:), pointer :: cpro_met_qv

! Initialisation

temp = t0
dens = ro0
press = dens*rair*temp ! ideal gas law
hspec = 0.0d0
met_qv_id = -1

if (imeteo.ge.2) then
  call field_get_val_s_by_name('meteo_temperature', cpro_met_temp)
  call field_get_val_s_by_name('meteo_pressure', cpro_met_p)
  call field_get_id_try('meteo_humidity', met_qv_id)
  if (met_qv_id.ge.0) then
    call field_get_val_s(met_qv_id, cpro_met_qv)
  endif
endif

if (photolysis) then
  iphotolysis = 1
else
  iphotolysis = 2
endif

if (ippmod(iatmos).ge.1) then
  call field_get_val_s(icrom, crom)
  call field_get_val_s(itempc, cpro_tempc)
endif

if (ippmod(iatmos).ge.2) then
  call field_get_val_s(ivarfl(isca(iymw)), cvar_totwt)
  call field_get_val_s(iliqwt, cpro_liqwt)
endif

! Computation of kinetic rates in every cell

! Computation of zenith angle
qureel = float(squant)
heurtu = float(shour) + float(smin)/60.d0+ssec/3600.d0
if (idtvar.eq.0 .or. idtvar.eq.1) heurtu = heurtu + ttcabs/3600.d0
call raysze(xlat,xlon,qureel,heurtu,0,albe,dlmuzero, omega, fo)
azi = dabs(dacos(dlmuzero)*180.d0/pi)

! Note: Photolysis should be cut in SPACK and not here even if the azimuthal angle is > 90

! Loop on cells
do iel = 1, ncel
  zent = xyzcen(3,iel) ! Z coordinate of the cell

  ! Temperature and density
  ! Dry or humid atmosphere
  if (ippmod(iatmos).ge.1) then
    temp = cpro_tempc(iel) + tkelvi
    dens = crom(iel)
    press = dens*rair*temp

    ! Constant density with a meteo file
  else if (imeteo.eq.1) then

    ! Hydrostatic pressure
    call intprf                                                   &
   (nbmett, nbmetm,                                               &
    ztmet , tmmet, phmet, zent, ttcabs, press )

    ! Temperature
    call intprf                                                   &
   (nbmett, nbmetm,                                               &
    ztmet , tmmet, ttmet, zent, ttcabs, temp )

    temp = temp + tkelvi
  else
    press = cpro_met_p(iel)
    temp = cpro_met_temp(iel) + tkelvi
  endif

  ! Specific humidity
  ! Humid atmosphere
  if (ippmod(iatmos).ge.2) then
    hspec = (cvar_totwt(iel)-cpro_liqwt(iel))    &
          / (1.d0-cpro_liqwt(iel))

  ! Constant density or dry atmosphere with a meteo file
  else if (imeteo.eq.1) then

    call intprf                                                   &
   (nbmett, nbmetm,                                               &
    ztmet , tmmet, qvmet, zent, ttcabs, hspec )

  else if (met_qv_id.ge.0) then
    hspec = cpro_met_qv(iel)

  ! Dry
  else
    hspec = 0.d0
  endif

  ! Call the computation of kinetic rates
  if (ichemistry.eq.1) then
    call kinetic_1(nrg,rk,temp,hspec,press,azi,1.0d0,iphotolysis)
  else if (ichemistry.eq.2) then
    call kinetic_2(nrg,rk,temp,hspec,press,azi,1.0d0,iphotolysis)
  else if (ichemistry.eq.3) then
    call kinetic_3(nrg,rk,temp,hspec,press,azi,1.0d0,iphotolysis)
  else if (ichemistry.eq.4) then
    call kinetic_4(nrg,rk,temp,hspec,press,azi,1.0d0,iphotolysis)
  endif

  ! Storage of kinetic rates
  do ii = 1, nrg
    reacnum((ii-1)*ncel+iel) = rk(ii)
  enddo

enddo

!--------
! Formats
!--------
!----
! End
!----

return
end subroutine kinrates
