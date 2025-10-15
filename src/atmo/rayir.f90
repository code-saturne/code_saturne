!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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
!> \brief Compute infrared flux divergence profile and downward flux at
!> ground level relying on a 1D radiative scheme.
!>
!> More precisely, compute atmospheric infrared (IR) radiation model quantities:
!> - vertical profile of IR flux divergence,
!> - downward IR flux at ground level
!> - upward and downward fluxes at different vertical levels
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ivertc      index of vertical profile
!> \param[in]     k1          index corresponding to ground level
!> \param[in]     kmray       number of vertical levels for radiation computation
!> \param[in]     emis        ground surface emissivity
!> \param[in,out] qqv         water vapor + dimers optical depth (0,zqq)
!> \param[in,out] qqqv        idem qqv but for intermediates vertical levels (zray)
!> \param[in]     qqvinf      idem qqv but for contribution above 11000m
!> \param[in]     zqq         vertical coordinate
!> \param[in]     zray        altitude (physical mesh)
!> \param[in]     temray      temperature in Celsius
!> \param[in]     qvray       specific humidity for water vapor
!> \param[in]     qlray       specific humidity for liquid water
!> \param[in]     fnerir      cloud fraction
!> \param[in]     romray      air density
!> \param[in]     preray      pressure
!> \param[in]     aeroso      aerosol concentration in micro-g/m3
!> \param[in]     t_surf      surface temperature
!> \param[in]     qw_surf     surface total water mass fraction
!> \param[in]     rho_surf    surface density
!> \param[in]     p_surf      surface pressure
!> \param[out]    foir        downward IR flux at the ground
!> \param[out]    rayi        IR flux divergence
!> \param[in]     ncray       Number of droplets interpolated on vertical grid
!-------------------------------------------------------------------------------

subroutine rayir &
 ( ivertc, k1, kmray, emis,                       &
   qqv, qqqv, qqvinf, zqq,                        &
   zray, temray, qvray, qlray, fnerir,            &
   romray, preray, aeroso,                        &
   t_surf, qw_surf, rho_surf, p_surf, foir, rayi, ncray )

!===============================================================================
! Module files
!===============================================================================

use optcal
use cstphy
use parall
use ppincl
use cs_c_bindings
use mesh
use field
use atincl, only: kmx, sigc, iru, ird, cp_a, cp_v, rad_atmo_model,zaero,aod_ir
use atincl, only: conco2
use cstnum, only: epzero, pi

!===============================================================================

implicit none

procedure() :: rayigc, rayive

! Arguments

integer ivertc, k1, kmray
double precision emis
double precision qqv(kmx+1), qqqv(kmx+1), qqvinf, zqq(kmx+1)
double precision temray(kmx), romray(kmx), preray(kmx)
double precision qvray(kmx), qlray(kmx), fnerir(kmx)
double precision zray(kmx)
double precision ncray(kmx)
double precision aeroso(kmx)
double precision t_surf, qw_surf, rho_surf, p_surf
double precision foir, rayi(kmx)

! Local variables

integer k, kk, ineb, inua, iaer, i, iaero_top
integer          ifac, iz1, iz2, f_id, c_id, iel
double precision sig
double precision qqcinf, cetyp, dzs8
double precision xqqvinf, xqqcinf, xqqlinf, fo, t4zt
double precision tauv, dtauv, taul, ctray, a1
double precision qqqqv, qqqqc, qqqql, a2, tvinfe
double precision dtvinfe, tvsup, dtvsup, dul
double precision t41, tlinfe, tlsup
double precision fn, fns, fni, u, tco2, zz, zzk
double precision pzp0, zqq0, corp
double precision beta, wh2ol
double precision zbas
double precision alpha_k
double precision rm, req
double precision a3, tvsups, dtvsups
double precision foirs, foirs1, foirs2
double precision tlsups, fnss
double precision var, zent, cpvcpa
double precision cort,xqqco2inf,abco2,dabco2,qqqqco2
double precision dz_aero
double precision abcsup,dabcsup,abcsups,dabcsups,abcinfe,dabcinfe
logical          is_active

integer, dimension(:), pointer :: itypfb
type(c_ptr) :: c_itypfb

double precision, allocatable :: rov(:), roc(:), rol(:), qv0(:), qc(:)
double precision, allocatable :: roco2(:),qco2(:),qqco2(:),qqqco2(:)
double precision, allocatable :: qqc(:), qql(:)
double precision, allocatable :: qqqc(:), qqql(:)
double precision, allocatable :: pspo(:), pspoqq(:), dt4(:), tqq(:)
double precision, allocatable :: dz0(:)
double precision, allocatable :: kliq(:)

double precision, allocatable :: dfir(:), ufir(:)
double precision, allocatable :: ckup(:), ckdown(:)
double precision, dimension(:,:), pointer :: bpro_rad_inc
double precision, dimension(:,:), pointer :: cpro_ck_up
double precision, dimension(:,:), pointer :: cpro_ck_down

!===============================================================================

allocate(rov(kmx), roc(kmx), rol(kmx), qv0(kmx), qc(kmx))
allocate(roco2(kmx), qco2(kmx))
allocate(qqc(kmx+1), qql(kmx+1))
allocate(qqco2(kmx+1))
allocate(qqqc(kmx), qqql(kmx))
allocate(qqqco2(kmx))
allocate(pspo(kmx), pspoqq(kmx), dt4(kmx), tqq(kmx))
allocate(dz0(kmx))
allocate(kliq(kmx+1))

! flux divided by Stefan constant for accuracy
allocate(dfir(kmx+1), ufir(kmx+1))
allocate(ckup(kmx), ckdown(kmx))

! local initializations

! indexes for presence of clouds, aerosols
inua = 0
iaer = 0
iaero_top = 0
dz_aero = 0.d0

cpvcpa = cp_v / cp_a
xqqlinf = 0.d0
taul = 0.d0
sig = stephn ! Boltzmann constant 5.6703d-8

do k = 1, kmx
  rov(k) = 0.d0
  roc(k) = 0.d0
  rol(k) = 0.d0
  roco2(k)=0.d0
  qv0(k) = 0.d0
  qc (k) = 0.d0
  qco2(k)= 0.d0
  qqv(k) = 0.d0
  qqc(k) = 0.d0
  qqco2(k)=0.d0
  qql(k) = 0.d0
  qqqv(k) = 0.d0
  qqqc(k) = 0.d0
  qqqco2(k)=0.d0
  qqql(k) = 0.d0
  pspo(k) = 0.d0
  pspoqq(k) = 0.d0
  dt4(k) = 0.d0
  tqq(k) = 0.d0
  rayi(k) = 0.d0

  if (aeroso(k).gt.1.d-8) iaer = 1
enddo

qqv(kmx+1) = 0.d0
qqc(kmx+1) = 0.d0
qqco2(kmx+1) = 0.d0
qql(kmx+1) = 0.d0
qqqv(kmx+1) = 0.d0

foir = 0.d0
! upper layers contribution (11000-44000) to optical depth
qqvinf = 0.005d0
qqcinf = 3.28d-7
cetyp = rvsra / 1.134d0

! diffusion coefficient integrated over all directions
beta = 1.66d0

! calculation for the liquid water absorption coefficient
! (same equivalent radius as solar radiation)

do k = k1, kmray
  ! liquid water content in g/m3 in the layers
  wh2ol = 1.d3*romray(k)*qlray(k)

  ! absorption coefficient depending on the equivalent radius
  ! of the cloud droplets
  rm = 30.d0 * wh2ol + 2.d0

  ! in case, microphysics data is available
  if (ncray(k).gt.epzero.and.qlray(k).gt.epzero) then
    req = 1.d6*( (3.d0*romray(k)*qlray(k)) /                     &
                 (4.d0*pi*1000.d0*ncray(k)*1.d6) )**(1./3.)      &
          *exp(sigc**2)
  ! otherwise
  else
    req = 1.5d0 * rm
  endif

  kliq(k) = beta * 0.75d0 * 1.d3 / req

  ! estimation of the inverse of the aerosol layer
  if((iaer.eq.1).and.(zray(k).le.zaero)) then
    iaero_top = max(k+1, iaero_top)
    fnerir(k) = 1.d0
  endif
enddo

if (iaero_top.ne.0) then
  dz_aero = 1.d0/zqq(iaero_top)
endif

do k = k1, kmray-1
  dz0(k) = zray(k+1) - zray(k)
enddo

dz0(kmray) = 0.d0
zbas = zqq(k1)
! Warning zqq(kmx+1) is not 16 000 for the moment...
zqq0 = 44000.d0

! 1. Computation to estimate absorption for water vapor and its dimer

do k = k1, kmray
  if(qlray(k).gt.1.d-8) inua = 1

  ! Note, simplification:
  ! so pspo(k) * preray(k1) = preray(k)
  corp = preray(k) / 101315.d0
  cort= tkelvi/(temray(k) + tkelvi)
  pspo(k) = preray(k) / p_surf

  qv0(k) = qvray(k)*corp*sqrt(cort)
  rov(k) = romray(k)*qv0(k)
  qco2(k)= conco2*(corp)**0.75d0*(cort)**0.325d0
  roco2(k)=romray(k)*qco2(k)
  qc(k) =   qvray(k)*qvray(k)*corp                               &
          * exp(1800.d0*(1.d0/(temray(k) + tkelvi)-1.d0/296.d0)) &
          * cetyp
  roc(k) = romray(k)*qc(k)
enddo

! idem for cloud layers (aerosol only for clear sky condition)
if(inua.eq.1) then

  do k = k1, kmray
    if (zray(k).le.zaero) then
      rol(k) = kliq(k)*romray(k)*qlray(k) + beta * aod_ir * dz_aero
    else
      rol(k) = kliq(k)*romray(k)*qlray(k)
    endif
  enddo

! idem for aerosol layers
elseif(iaer.eq.1) then

  do k = k1, kmray
    if (zray(k).le.zaero) then
      rol(k) = beta * aod_ir * dz_aero
    endif
  enddo

endif

! 2. optical depth calculation for water vapor and its dimer
! qqq corresponds to standard levels, qq to intermediate levels

! surface temperature (Kelvin)
tqq(k1) = t_surf + tkelvi

do k = k1+1, kmray
  alpha_k = ( zray(k) - zqq(k)) / (zray(k) - zray(k-1))
  tqq(k) = (alpha_k * temray(k)+(1.d0 - alpha_k)*temray(k-1)) + tkelvi
  pspoqq(k) = (alpha_k * pspo(k) + (1.d0 - alpha_k) * pspo(k-1))

  ! dT^4/dz x dz = 4 T^3 dT/dz x dz
  dt4(k) = 4.d0 * (tqq(k)**3) * (temray(k) - temray(k-1))
enddo

! qqv(k+1) = int(0, zqq(k+1)) rov dz
! qqqv(k) = qqv(k) + dz/2(k) * rov(k)
do k = k1, kmray
  qqqv(k) = qqv(k) + (zray(k) - zqq(k)) * rov(k)
  qqqc(k) = qqc(k) + (zray(k) - zqq(k)) * roc(k)
  ! int(0, zray(k))
  qqqco2(k) = qqco2(k) + (zray(k) - zqq(k)) * roco2(k)
  qqv(k+1) = qqqv(k) + (zqq(k+1) - zray(k)) * rov(k)
  qqc(k+1) = qqqc(k) + (zqq(k+1) - zray(k)) * roc(k)
  ! int(0, zqq(k+1))
  qqco2(k+1) = qqqco2(k) + (zqq(k+1) - zray(k)) * roco2(k)
enddo

pspoqq(k1) = pspo(k1)
! TODO should be qqv(kmray+1)
xqqvinf = qqqv(kmray) + qqvinf
xqqcinf = qqqc(kmray) + qqcinf
! Global integral from 0 to 4km
xqqco2inf = qqco2(kmray+1)

! 3. optical depth calculation for liquid water in cloudy cases

if(inua.eq.1) then
  do k = k1, kmray
    ! int(0, zray(k))
    qqql(k) = qql(k) + (zray(k) - zqq(k)) * rol(k)
    ! int(0, zqq(k+1))
    qql(k+1) = qqql(k) + (zqq(k+1) - zray(k)) * rol(k)
  enddo

  xqqlinf = qql(kmray+1)
endif

! 5. IR downward flux computation for the ground level

fo = 0.d0
t4zt = (temray(kmray) + tkelvi)**4
t41 = (tqq(k1))**4

! For cloudy sky (but also valid for clear sky)
fn = fnerir(k1+1)
do k = k1+1, kmray
  ! cloud fraction estimation (we take the maximum)
  fn = max(fn,fnerir(k))
  if (fn.lt.1.d-3) then
    fn = 0.d0
    taul = 1.d0
  else
    taul = exp(-qql(k)/fn)
  endif

  call rayive(tauv,dtauv,qqv(k),qv0(k1),qqc(k),qc(k1),romray(k))
  ! Compute absorption from zbas to zqq(k)
  call rayigc(zqq(k1), zqq(k), abco2,dabco2,qqv(k),qv0(k1),qqco2(k),qco2(k1),romray(k))
  fo = fo - (1.d0 - (1.d0 + fn*(taul - 1.d0))                     &
          * (tauv - abco2)) * dt4(k)
enddo

! Last level: 11km to 44km, assumed to be isothermal
call rayive(tauv,dtauv,xqqvinf,qv0(k1),xqqcinf,qc(k1),            &
            romray(kmray))

call rayigc(zqq(k1), zqq0, abco2,dabco2,xqqvinf,qv0(k1),xqqco2inf,qco2(k1), &
            romray(kmray))
fo = fo + t4zt*(1.d0 - (1.d0 + fn*(taul - 1.d0))*(tauv - abco2))

foir = fo
! IR flux calculation in the vertical layers
do i = k1, kmray
  dfir(i) = 0.d0
  ufir(i) = 0.d0
  fn = fnerir(i)
  do k = i+1, kmray
    ! cloud fraction estimation (we take the maximum)
    fn = max(fn,fnerir(k))
    qqqqv = qqv(k)-qqqv(i)
    qqqqc = qqc(k)-qqqc(i)
    qqqqco2 = qqco2(k)-qqqco2(i)
    qqqql = qql(k)-qqql(i)

    if (fn.lt.1.d-3) then
      fn = 0.d0
      taul = 1.d0
    else
      taul = exp(-qqqql/fn)
    endif

    ! downward fluxes
    call rayive(tauv, dtauv, qqqqv, qv0(i), qqqqc, qc(i), romray(i))
    ! Compute from i to level k
    call rayigc(zqq(i), zqq(k), abco2,dabco2,qqqqv,qv0(i),qqqqco2,qco2(i),romray(i))
    dfir(i) = dfir(i)-(1.d0-(1.d0+fn*(taul-1.d0))*(tauv-abco2)) * dt4(k)
  enddo

  qqqqv = xqqvinf-qqqv(i)
  qqqqc = xqqcinf-qqqc(i)
  qqqqco2 = xqqco2inf-qqqco2(i)
  call rayive(tauv, dtauv, qqqqv, qv0(i), qqqqc, qc(i), romray(i))
  call rayigc(zqq(k1), zqq0, abco2,dabco2,qqqqv,qv0(i),qqqqco2,qco2(i),romray(i))

  dfir(i) = dfir(i)+t4zt*(1.d0-(1.d0+fn*(taul-1.d0))*(tauv-abco2))

  ! upward fluxes
  if (i.gt.k1) then
    fn = fnerir(i)
    do k = k1+1, i
      ! cloud fraction estimation (we take the maximum)
      fn = max(fn,fnerir(k))
      qqqqv = qqqv(i)-qqv(k)
      qqqqc = qqqc(i)-qqc(k)
      qqqqco2 = qqqco2(i)-qqco2(k)
      qqqql = qqql(i)-qql(k)

      if (fn.lt.1.d-3) then
        fn = 0.d0
        taul = 1.d0
      else
        taul = exp(-qqqql/fn)
      endif

      call rayive(tauv, dtauv, qqqqv, qv0(i), qqqqc, qc(i), romray(i))
      call rayigc(zqq(i), zqq(k), abco2,dabco2,qqqqv,qv0(i),qqqqco2,qco2(i),romray(i))

      ufir(i) = ufir(i)+(1.d0-(1.d0+fn*(taul-1.d0))*(tauv-abco2))*dt4(k)
    enddo
  endif

  ! contribution of the upward flux reflected at the ground
  ! The modification proposed by Ponnulakshmi and al. (2009) for the upward
  ! flux reflected by the ground when emissivity is different from 1
  ! requires to compute the integral below (a3).
  a3 = 0.d0
  do k = k1+1, kmray
    qqqqv = qqqv(i) + qqv(k)
    qqqqc = qqqc(i) + qqc(k)
    qqqqco2 = qqqco2(i) + qqco2(k)
    qqqql = qqql(i) + qql(k)

    call rayive(tauv, dtauv, qqqqv, qv0(i), qqqqc, qc(i), romray(i))
    ! Integration from i to k
    call rayigc(zqq(i), zqq(k), abco2,dabco2,qqqqv,qv0(i),qqqqco2,qco2(i),romray(i))
    if(fn.lt.1.d-3) then
      fn = 0.d0
      taul = 1.d0
    else
      taul = exp(-qqqql/fn)
    endif
    a3 = a3 - (1.d0-(1.d0+fn*(taul-1.d0))*(tauv-abco2)) * dt4(k)
  enddo

  qqqqv = qqqv(i) + xqqvinf
  qqqqc = qqqc(i) + xqqcinf
  qqqqco2 = qqqco2(i) + xqqco2inf

  call rayive(tvsups, dtvsups, qqqqv, qv0(i), qqqqc, qc(i), romray(i))
  ! Integration from k1 to i
  call rayigc(zqq(k1),zqq(i), abcsups,dabcsups,qqqqv,qv0(i),qqqqco2,qco2(i),romray(i))

  foirs = a3+(1.d0-(1.d0+fn*(taul-1.d0))*(tvsups-abcsups))*t4zt

  ! Upward fluxes estimation (sum of direct part and reflected part)
  if (i.gt.k1) then
    ufir(i) = ufir(i) + emis*t41 + (1.d0-emis)*foirs
  else
    ufir(k1) = emis*t41 + (1.d0-emis)*foir
  endif
enddo

! 6. Cooling in the vertical layers

! For cloudy sky (also valid for clear sky)
do k = k1, kmray

  ctray = sig/romray(k)/(cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(k)))
  ! a1: contribution from 0 to z
  a1 = 0.d0
  dul = rol(k)

  do kk = k1+1, k
    qqqqv = qqqv(k) - qqv(kk)
    qqqqc = qqqc(k) - qqc(kk)
    qqqqco2 = qqqco2(k) - qqco2(kk)
    qqqql = qqql(k) - qql(kk)

    ! Integration from k to kk
    call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    call rayigc(zqq(k), zqq(kk), abco2,dabco2,qqqqv,qv0(k),qqqqco2,qco2(k),romray(k))

    fn = fnerir(kk)

    do ineb = kk, k
      fn = max(fn,fnerir(ineb))
    enddo

    if(fn.lt.1.d-3) then
      fn = 0.d0
      taul = 1.d0
    else
      taul = exp(-qqqql/fn)
    endif

    a1 = a1 + dt4(kk)*(  (1.d0 + fn*(taul - 1.d0))*(dabco2 - dtauv)   &
                         + taul*dul*(tauv - abco2))
  enddo

  ! a2: contribution from z to ztop
  a2 = 0.d0
  do kk = k+1, kmray
    qqqqv = qqv(kk) - qqqv(k)
    qqqqc = qqc(kk) - qqqc(k)
    qqqqco2 = qqco2(kk) -  qqqco2(k)
    qqqql = qql(kk) - qqql(k)

    call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    ! Integration from k to kk
    call rayigc(zqq(k), zqq(kk), abco2,dabco2,qqqqv,qv0(k),qqqqco2,qco2(k),romray(k))

    fn = fnerir(k)

    do ineb = k, kk
      fn = max(fn,fnerir(ineb))
    enddo

    if (fn.lt.1.d-3) then
      fn = 0.d0
      taul = 1.d0
    else
      taul = exp(-qqqql/fn)
    endif
    a2 = a2 + dt4(kk)*(  (1.d0 + fn*(taul - 1.d0))*(dabco2-dtauv)    &
                         + taul*dul*(tauv - abco2))
  enddo

  ! a3: contribution of the upward flux reflected at the ground
  ! Modification of Ponnulakshmi and al. (2009)
  a3 = 0.d0
  do kk = k1+1, kmray
    qqqqv = qqqv(k) + qqv(kk)
    qqqqc = qqqc(k) + qqc(kk)
    qqqqco2 = qqqco2(k) + qqco2(kk)
    qqqql = qqql(k) + qql(kk)

    call rayive(tauv, dtauv, qqqqv, qv0(k), qqqqc, qc(k), romray(k))
    ! Integration from k to kk
    call rayigc(zqq(k), zqq(kk), abco2,dabco2,qqqqv,qv0(k),qqqqco2,qco2(k),romray(k))
    do ineb = k1, kmray
      fn = max(fn,fnerir(ineb))
    enddo

    if(fn.lt.1.d-3) then
      fn=0.d0
      taul=1.d0
    else
      taul=exp(-qqqql/fn)
    endif

    a3 = a3 + dt4(kk)*( (1.d0+fn*(taul-1.d0))*(dabco2-dtauv)           &
                         +taul*dul*(tauv-abco2))
  enddo

  ! contribution from z to infinity
  qqqqv = xqqvinf - qqqv(k)
  qqqqc = xqqcinf - qqqc(k)
  qqqqco2 = xqqco2inf - qqqco2(k)
  qqqql = xqqlinf - qqql(k)
  call rayive(tvsup, dtvsup, qqqqv, qv0(k), qqqqc, qc(k), romray(k))
  call rayigc(zqq(k), zqq0, abcsup,dabcsup,qqqqv,qv0(k),qqqqco2,qco2(k),romray(k))

  fns = fnerir(k)
  do ineb = k, kmray
    fns = max(fns,fnerir(ineb))
  enddo

  if(fns.lt.1.d-3) then
    fns = 0.d0
    tlsup = 1.d0
  else
    tlsup = exp(-qqqql/fns)
  endif

  !  ground contribution transmitted by lower layers (0-z)
  qqqqv = qqqv(k)
  qqqqc = qqqc(k)
  qqqqco2 = qqqco2(k)
  qqqql = qqql(k)
  call rayive(tvinfe, dtvinfe, qqqqv, qv0(k), qqqqc, qc(k), romray(k))
  ! Integration from
  call rayigc(zqq(k1), zqq(k),abcinfe,dabcinfe,qqqqv,qv0(k),qqqqco2,qco2(k),romray(k))

  fni = fnerir(k1+1)
  do ineb =2, k
    fni = max(fni,fnerir(ineb))
  enddo

  if(fni.lt.1.d-3) then
    tlinfe = 0.d0
    fni = 0.d0
  else
    tlinfe=exp(-qqqql/fni)
  endif
  ! contribution of the upward flux reflected by the ground
  qqqqv = qqqv(k)+xqqvinf
  qqqqc = qqqc(k)+xqqcinf
  qqqqco2 = qqqco2(k)+xqqco2inf
  qqqql = qqql(k)+xqqlinf
  call rayive(tvsups, dtvsups, qqqqv, qv0(k), qqqqc, qc(k), romray(k))
  call rayigc(zqq(k1), zqq0, abcsups,dabcsups,qqqqv,qv0(k),qqqqco2,qco2(k),romray(k))

  fnss=fnerir(k)
  do ineb=k1,kmray
    fnss=max(fnss,fnerir(ineb))
  enddo

  if(fnss.lt.1.d-3) then
    tlsups=0.d0
    fnss=0.d0
  else
    tlsups=exp(-qqqql/fnss)
  endif

  ! cooling rate in cloudy conditions
  ! Formula follows Ponnulakshmi and al. (2009)
  ! TODO should give back previous formula for dul=0 and tlsup=1 and tlsups=1
  rayi(k) =  ctray*(a1 - a2 + t4zt*((1.d0 + fns*(tlsup - 1.d0))             &
           *(dabcsup - dtvsup) + dul*tlsup*(tvsup - abcsup))                &
           -(1.d0-emis)*(a3 + t4zt*((1.d0 + fnss*(tlsups-1.d0))             &
           *(dtvsups-dabcsups)-dul*tlsups*(tvsups-abcsups))))

  ! Save for 3D
  ! Ck_up = depsg0 / (1 - epsg0)
  ! With:
  ! epsg0 = 1 - tv_infinity + acinfinity
  ! depsg0 = - dtv_infinity + dacinfinity
  ckup(k) = (dabcinfe - dtvinfe) / (tvinfe - abcinfe)

  ! Ck_down = -depsgz / (1 - epsgz)
  ! With:
  ! epsgz = 1 - tv_sup + acinfinity
  ! depsgz = - dtv_infinity + dacinfinity
  ckdown(k) = (dabcsup - dtvsup) / (tvsup - abcsup)

enddo

! Finalization: multiplication by sig
foir = sig * foir
do k = 1, kmray
  ufir(k) = sig * ufir(k)
  dfir(k) = sig * dfir(k)
  iru(k,ivertc) = ufir(k)
  ird(k,ivertc) = dfir(k)
enddo

! Compute Boundary conditions for the 3D InfraRed radiance
! at the top of the CFD domain
call field_get_id_try("spectral_rad_incident_flux", f_id)

is_active = cs_rad_time_is_active()

if (f_id.ge.0.and.is_active) then
  call field_get_val_v(f_id, bpro_rad_inc)

  call field_get_val_v_by_name("rad_absorption_coeff_up", cpro_ck_up)
  call field_get_val_v_by_name("rad_absorption_coeff_down", cpro_ck_down)

  ! First count if solar is activated
  c_id = 0
  ! Direct solar radiation incident (for H2O band)
  if (iand(rad_atmo_model, 1).eq.1) then
    c_id = c_id + 1
  endif

  ! Direct solar radiation incident (for O3 band)
  if (iand(rad_atmo_model, 2).eq.2) then
    c_id = c_id + 1
  endif

  ! Diffuse solar radiation incident
  if (iand(rad_atmo_model, 4).eq.4) then
    c_id = c_id + 1
  endif

  ! Diffuse solar radiation incident (for SUV O3 band)
  if (iand(rad_atmo_model, 8).eq.8) then
    c_id = c_id + 1
  endif

  ! Infra Red radiation incident
  if (iand(rad_atmo_model, 16).eq.16) then

    call cs_f_boundary_conditions_get_pointers(c_itypfb)
    call c_f_pointer(c_itypfb, itypfb, [nfabor])

    c_id = c_id + 1
    do ifac = 1, nfabor

      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else
        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          dfir, zent, iz1, iz2, var )

        bpro_rad_inc(c_id, ifac) = var
      endif
    enddo

    ! Store the (downward and upward) absorption coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      call intprz &
        (kmray, zray,                                               &
        ckdown, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = var * 3.d0 / 5.d0

      call intprz &
        (kmray, zray,                                               &
        ckup, zent, iz1, iz2, var )

      cpro_ck_up(c_id, iel) = var * 3.d0 / 5.d0

    enddo

  endif

endif

deallocate(ufir,dfir)
deallocate(ckup, ckdown)

deallocate(rov,roc,rol,qv0,qc)
deallocate(roco2, qco2)
deallocate(qqc,qql)
deallocate(qqco2)
deallocate(qqqc,qqql)
deallocate(qqqco2)
deallocate(pspo,pspoqq,dt4,tqq)
deallocate(dz0)
deallocate(kliq)

return
end subroutine rayir
