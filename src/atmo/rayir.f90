!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \file rayir.f90
!> \brief 1D Radiative scheme - IR flux divergence profile and
!>  downward flow to the ground
!
!>  \brief   calculation for Infrared (IR)  atmospheric radiation
!>      - vertical profile of IR flux divergence,
!>      - downward IR flux at the ground.
!>      - upward and downward fluxes at different vertical levels for irdm=1 (option)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]   k1          index corresponding to ground level
!> \param[in]   kmray       number of vertical levels for radiation computation
!> \param[in]   ico2        ico2=1 -> compute CO2 absorption
!> \param[in]   emis        ground surface emissivity
!> \param[in]   qqv         water vapor + dimers optical depth (0,z)
!> \param[in]   qqqv        idem qqv but for intermediates vertical levels
!> \param[in]   qqvinf      idem qqv but for contribution above 11000m
!> \param[in]   zqq         vertical coordinate
!> \param[in]   acinfe      absorption for CO2 + O3 (0,z)
!> \param[in]   dacinfe     differential absorption for CO2 + 03 (0,z)
!> \param[in]   aco2        idem acinfe but for CO2 only
!> \param[in]   daco2       idem dacinfe but for CO2 only
!> \param[in]   acsup       idem acinfe, for (z,0)
!> \param[in]   dacsup      idem dacinfe, for (z,0)
!> \param[in]   zray        altitude (physical mesh)
!> \param[in]   temray      temperature in Celsius
!> \param[in]   qvray       specific humidity for water vapor
!> \param[in]   qlray       specific humidity for liquid water
!> \param[in]   fnerir      cloud fraction
!> \param[in]   romray      air density
!> \param[in]   preray      pressure
!> \param[in]   aeroso      aerosol concentration in micro-g/m3
!> \param[out]  foir        downward IR flux at the ground
!> \param[out]  rayi        IR flux divergence
!-------------------------------------------------------------------------------

subroutine rayir &
 ( k1,kmray,ico2,emis,                            &
   qqv, qqqv, qqvinf, zqq,                        &
   acinfe, dacinfe, aco2, daco2, acsup, dacsup,   &
   zray,temray,qvray,qlray,fnerir,                &
   romray,preray,aeroso,foir,rayi )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl, only: kmx, rvsra, cpvcpa

!===============================================================================

implicit none

! Arguments

integer k1,kmray,ico2
double precision emis
double precision qqv(kmx+1), qqqv(kmx+1), qqvinf, zqq(kmx)
double precision acinfe(kmx), dacinfe(kmx), aco2(kmx,kmx), daco2(kmx,kmx)
double precision acsup(kmx), dacsup(kmx)
double precision temray(kmx),romray(kmx),preray(kmx)
double precision qvray(kmx),qlray(kmx),fnerir(kmx)
double precision zray(kmx)
double precision aeroso(kmx)
double precision foir,rayi(kmx)

!Local variables

integer n,k,kk,kp1,ineb,inua,iaer
double precision sig
double precision qqcinf,cetyp,dzs8
double precision xqqvinf,xqqcinf,xqqlinf,fo,t4zt
double precision tauv,dtauv,taul,ctray,a1
double precision qqqqv,qqqqc,qqqql,a2,tvinfe
double precision dtvinfe,tvsup,dtvsup,dul
double precision t41,tlinfe,tlsup,qqlinf
double precision fn,fns,fni,u,tco2,zz,zzk
double precision pzp0,zqq0,corp
double precision beta,wh2ol
double precision zbas
double precision caero

double precision,allocatable :: rov(:),roc(:),rol(:),qv0(:),qc(:)
double precision,allocatable :: qqc(:),qql(:)
double precision,allocatable :: qqqc(:),qqql(:)
double precision,allocatable :: pspo(:),pspoqq(:),dt4dz(:),tqq(:)
double precision,allocatable :: dz0(:)
double precision,allocatable :: kliq(:)

!===============================================================================

allocate(rov(kmx),roc(kmx),rol(kmx),qv0(kmx),qc(kmx))
allocate(qqc(kmx),qql(kmx))
allocate(qqqc(kmx),qqql(kmx))
allocate(pspo(kmx),pspoqq(kmx),dt4dz(kmx),tqq(kmx))
allocate(dz0(kmx))
allocate(kliq(kmx))

!   0- initialisations locales
!   --------------------------
xqqlinf = 0.d0
taul = 0.d0
sig = stephn ! Bolztman constant 5.6703d-8

do n = 1, kmx
  rov(n) = 0.d0
  roc(n) = 0.d0
  rol(n) = 0.d0
  qv0(n) = 0.d0
  qc (n) = 0.d0
  qqv(n) = 0.d0
  qqc(n) = 0.d0
  qql(n) = 0.d0
  qqqv(n) = 0.d0
  qqqc(n) = 0.d0
  qqql(n) = 0.d0
  pspo(n) = 0.d0
  pspoqq(n) = 0.d0
  dt4dz(n) = 0.d0
  tqq(n) = 0.d0
  rayi(n) = 0.d0
enddo
qqv(kmx+1) = 0.d0
qqqv(kmx+1) = 0.d0

foir = 0.d0

! uppper layers contribution (11000-44000) to optical depth

qqvinf = 0.002515d0
qqcinf = 2.98d-7
qqlinf = 0.d0

cetyp = rvsra/1.134d0

! diffusion coeficient integrated over all directions

beta = 1.66d0

! calculation for the liquid water absorption coefficient
!  (same equivalent radius as solar radiation)
!
do k = k1, kmray
! liquid water content in g/m3 in the layers
  wh2ol = 1.d3*romray(k)*qlray(k)
! absorptiion coefficient depending on the equivalent radius
! of the cloud droplets
  kliq(k) = 120.d0
enddo
! indexes for presence of clouds, aerosols
inua = 0
iaer = 0
! index to compute or not the downward and upward IR fluxes
!irdm=1
! constant for aerosol concentration which has to be in µg/m3
caero=1.d-9

do k = k1, kmray-1
 dz0(k) = zray(k+1) - zray(k)
enddo

dz0(kmray) = 0.d0
zbas = zray(k1)

!   1- computation to estimate absorption for water vapor and its dimer
!   --------------------------------------------------------------------
do k = k1, kmray
  if(qlray(k).gt.1.d-8) inua = 1
  if(aeroso(k).gt.1.d-8) iaer = 1
  pspo(k) = preray(k)/preray(k1)
  corp = pspo(k)*preray(k1)/101300.d0
  qv0(k) = qvray(k)*corp*sqrt(tkelvi/(temray(k) + tkelvi))
  rov(k) = romray(k)*qv0(k)
  qc(k) = qvray(k)*qvray(k)*corp*exp(1800.d0*                         &
      (1.d0/(temray(k) + tkelvi)-1.d0/296.d0))*cetyp
  roc(k) = romray(k)*qc(k)
enddo

! idem for cloud layers (aerosol only for clear sky condition)
if(inua.eq.1) then
  do k = k1, kmray
    rol(k) = romray(k)*(qlray(k) + caero*aeroso(k))
  enddo

! idem for aerosol layers
elseif(iaer.eq.1) then
  do k = k1, kmray
    rol(k) = caero*romray(k)*aeroso(k)
    if(aeroso(k).gt.1.d-10) fnerir(k) = 1.d0
  enddo
endif

!  2 - optcal depth calculation for water vapor and its dimer
!  ------------------------------------------------------------------
! qqq corresponds to standard levels, qq to intermediate levels

do k = k1+1, kmray
  dzs8 = dz0(k-1)/8.d0
  tqq(k) = (temray(k)+temray(k-1))/2.d0 + tkelvi
  zqq(k) = (zray(k) + zray(k-1))/2.d0
  dt4dz(k) = 4.d0*(tqq(k)**3.d0)*(temray(k) - temray(k-1))/dz0(k-1)
  pspoqq(k) = (pspo(k) + pspo(k-1))/2.d0
  qqv(k) = qqqv(k-1) + dzs8*(rov(k) + 3.d0*rov(k-1))
  qqc(k) = qqqc(k-1) + dzs8*(roc(k) + 3.d0*roc(k-1))
  qqqv(k) = qqv(k) + dzs8*(3.d0*rov(k) + rov(k-1))
  qqqc(k) = qqc(k) + dzs8*(3.d0*roc(k) + roc(k-1))
enddo

zqq(k1) = zray(k1)
pspoqq(k1) = pspo(k1)
xqqvinf = qqqv(kmray) + qqvinf
xqqcinf = qqqc(kmray) + qqcinf

!   3 - optical depth calculation for liquid water in cloudy cases
!       by using the same notation
!   ----------------------------------------------------------------

if(inua.eq.1) then
  do k = k1+1,kmray
    dzs8 = dz0(k-1)/8.d0
    qql(k) = qqql(k-1) + dzs8*(rol(k) + 3.d0*rol(k-1))
    qqql(k) = qql(k) + dzs8*(3.d0*rol(k) + rol(k-1))
  enddo
  xqqlinf = qqql(kmray) + qqlinf
endif

!   4 - in order to economize computation time, CO2 and O3 absorption are calculated only
!       at the first time step and stocked in tables

if (ico2.eq.1) then

  do k = k1,kmray
    zz = zray(k)
    if(k.ne.k1) then
      do kk = k1+1, k
        u = qqqv(k) - qqv(kk)
        tco2 = (temray(k) + tkelvi + tqq(kk))/2.d0
        call rayigc(zbas,zz,pspo(k),zqq(kk),pspoqq(kk),aco2(k,kk),    &
               daco2(k,kk),qv0(k),u,tco2,romray(k))
      enddo
    endif
    kp1 = k+1
    do kk = kp1, kmray
      u = qqv(kk) - qqqv(k)
      tco2 = (temray(k) + tkelvi + tqq(kk))/2.d0
      call rayigc(zbas,zz,pspo(k),zqq(kk),pspoqq(kk),aco2(k,kk),      &
             daco2(k,kk),qv0(k),u,tco2,romray(k))
    enddo
    u = xqqvinf - qqqv(k)
    tco2 = 375.d0
    pzp0 = 0.d0
    zqq0 = 44000.d0
    call rayigc(zbas,zz,pspo(k),zqq0,pzp0,acsup(k),dacsup(k),         &
           qv0(k),u,tco2,romray(k))
    if(k.ne.k1) then
      u = qqqv(k)
      tco2 = (temray(k1)+temray(k))/2.d0 + tkelvi
      zz = zray(k1)
      zzk = zray(k)
      call rayigc(zbas,zz,pspo(k1),zzk,pspo(k),acinfe(k),dacinfe(k),  &
             qv0(k),u,tco2,romray(k))
    endif
  enddo

endif

!   5 - IR downward flux computation for the ground level:foir
!   ----------------------------------------------

fo = 0.d0
t4zt = (temray(kmray) + tkelvi)**4

if(inua.ne.1) then

!  for clear sky
!  =================
  do k = k1+1, kmray
! transmissivity for water vapor and its dimer
    call rayive(tauv,dtauv,qqv(k),qv0(k1),qqc(k),qc(k1),romray(k))
    fo = fo - (1.d0 - tauv + aco2(k1,k))*dt4dz(k)*dz0(k-1)
  enddo
  call rayive(tauv,dtauv,xqqvinf,qv0(k1),xqqcinf,qc(k1),            &
           romray(kmray))
  foir = sig*(fo + t4zt*(1.d0 - tauv + acsup(k1)))

else
!  for cloudy sky
!  ===================
  fn = fnerir(k1+1)
  do k = k1+1, kmray
! cloud fraction estimation (we take the maximum)
    fn = max(fn,fnerir(k))
    if(fn.lt.1.d-3) then
      fn = 0.d0
      taul = 0.d0
    else
      taul = exp(-kliq(k)*qql(k)/fn)
    endif
    call rayive(tauv,dtauv,qqv(k),qv0(k1),qqc(k),qc(k1),romray(k))
    fo = fo - (1.d0 - (1.d0 + fn*(taul - 1.d0))                                 &
       * (tauv - aco2(k1,k)))*dt4dz(k)*dz0(k-1)
  enddo
  call rayive(tauv,dtauv,xqqvinf,qv0(k1),xqqcinf,qc(k1),            &
           romray(kmray))
  foir = sig*(fo + t4zt*(1.d0 - (1.d0 + fn*(taul - 1.d0))*(tauv - acsup(k1))))
endif

!   6 -IR cooling in the vertcal layers
!   ----------------------------

t41 = (temray(k1) + tkelvi)**4
if(inua.ne.1) then

!  for clear sky
!  ==============
  do k = k1+1, kmray

    ctray = sig/romray(k)/(cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(k)))
    !  a1: contribution from  0 to z
    a1 = 0.d0
    do kk = k1+1, k
      qqqqv = qqqv(k) - qqv(kk)
      qqqqc = qqqc(k) - qqc(kk)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      a1 = a1 + dt4dz(kk)*(daco2(k,kk) - dtauv)*dz0(kk-1)
    enddo
    kp1 = k+1
    ! a2: contribution from z to ztop
    a2 = 0.d0
    do kk = kp1, kmray
      qqqqv = qqv(kk) - qqqv(k)
      qqqqc = qqc(kk) - qqqc(k)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      a2 = a2 + dt4dz(kk)*(daco2(k,kk) - dtauv)*dz0(kk-1)
    enddo
    ! contribution from z to infinity
    qqqqv = xqqvinf - qqqv(k)
    qqqqc = xqqcinf - qqqc(k)
    call rayive(tvsup,dtvsup,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    ! contribution du sol transmise par les couches inferieures (0-z)
    qqqqv = qqqv(k)
    qqqqc = qqqc(k)
    call rayive(tvinfe,dtvinfe,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    ! compute cooling
    rayi(k) = ctray*(a1 - a2 + t4zt*(dacsup(k) - dtvsup) -(1.d0 - emis)         &
              *(t41 - foir/sig)*(dtvinfe - dacinfe(k)))
  enddo
else

!  for cloudy sky
!  ================
  do k = k1+1, kmray

    ctray = sig/romray(k)/(cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(k)))
! a1: contribution from 0 to z
    a1 = 0.d0
    dul = kliq(k)*rol(k)
    do kk = k1+1, k
      qqqqv = qqqv(k) - qqv(kk)
      qqqqc = qqqc(k) - qqc(kk)
      qqqql = qqql(k) - qql(kk)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      fn = fnerir(kk)
      do ineb = kk, k
        fn = max(fn,fnerir(ineb))
      enddo
      if(fn.lt.1.d-3) then
        taul = 0.d0
        fn = 0.d0
      else
        taul = exp(-kliq(k)*qqqql/fn)
      endif
      a1 = a1 + dt4dz(kk)*((1.d0 + fn*(taul - 1.d0))*(daco2(k,kk) - dtauv)      &
         + taul*dul*(tauv - aco2(k,kk)))*dz0(kk-1)
    enddo
    kp1 = k+1
! a2: contribution from z to ztop
    a2 = 0.d0
    do kk = kp1, kmray
      qqqqv = qqv(kk) - qqqv(k)
      qqqqc = qqc(kk) - qqqc(k)
      qqqql = qql(kk) - qqql(k)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      fn = fnerir(k)
      do ineb = k, kk
        fn = max(fn,fnerir(ineb))
      enddo
      if(fn.lt.1.d-3) then
        taul = 0.d0
        fn = 0.d0
      else
        taul = exp(-kliq(k)*qqqql/fn)
      endif
      a2 = a2 + dt4dz(kk)*((1.d0 + fn*(taul - 1.d0))*(daco2(k,kk)-dtauv)        &
         + taul*dul*(tauv - aco2(k,kk)))*dz0(kk-1)
    enddo
! contribution de z a l'infini
    qqqqv = xqqvinf - qqqv(k)
    qqqqc = xqqcinf - qqqc(k)
    qqqql = xqqlinf - qqql(k)
    call rayive(tvsup,dtvsup,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    fns = fnerir(k)
    do ineb = k, kmray
      fns = max(fns,fnerir(ineb))
    enddo
    if(fns.lt.1.d-3) then
      tlsup = 0.d0
      fns = 0.d0
    else
      tlsup = exp(-kliq(k)*qqqql/fns)
    endif
!  ground contribution transmited by lower layers (0-z)
    qqqqv = qqqv(k)
    qqqqc = qqqc(k)
    qqqql = qqql(k)
    call rayive(tvinfe,dtvinfe,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    fni = fnerir(k1+1)
    do ineb =2, k
      fni = max(fni,fnerir(ineb))
    enddo
    if(fni.lt.1.d-3) then
      tlinfe = 0.d0
      fni = 0.d0
    else
      tlinfe=exp(-kliq(k)*qqqql/fni)
    endif
! calcul du refoidissement
    rayi(k) = ctray*(a1 - a2 + t4zt*((1.d0 + fns*(tlsup - 1.d0))                &
            *(dacsup(k) - dtvsup) + dul*tlsup*(tvsup - acsup(k)))               &
            -(1.d0 - emis)*(t41 - foir/sig)*((1.d0 + fni*(tlinfe - 1.d0))       &
            *(dtvinfe - dacinfe(k)) - dul*tlinfe*(tvinfe - acinfe(k))))

  enddo

endif

deallocate(rov,roc,rol,qv0,qc)
deallocate(qqc,qql)
deallocate(qqqc,qqql)
deallocate(pspo,pspoqq,dt4dz,tqq)
deallocate(dz0)
deallocate(kliq)

return
end subroutine rayir
