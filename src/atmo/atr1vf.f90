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
!> \file atr1vf.f90
!>
!> \brief Compute radiative fluxes for the atmospheric model.
!> Computes the source term for scalar equations from radiative forcing
!> (UV and IR radiative fluxes) with a 1D scheme.
!-------------------------------------------------------------------------------

subroutine atr1vf

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum, only:pi
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use atincl
use atsoil

!===============================================================================

implicit none

! Local variables
integer k, ii, jj
integer k1
integer ico2,imer1
integer ideb, icompt
integer kmray, ktamp

double precision heuray, albedo, emis, foir, fos
double precision xvert, yvert
double precision zrac,fpond,rap,tmoy,rhum,dum

integer , allocatable :: cressm(:), interp(:)
double precision, allocatable, dimension(:) :: temray, qvray, qlray, ncray
double precision, allocatable, dimension(:) :: fneray, romray, preray
double precision, allocatable, dimension(:) :: zproj, ttvert, qvvert, romvert
double precision, allocatable, dimension(:) :: qwvert, qlvert, ncvert,fnvert
double precision, allocatable, dimension(:) :: aeroso, infrad
double precision, allocatable, dimension(:,:,:) :: coords(:,:,:)
double precision, dimension(:), pointer :: crom, cpro_pcliq
double precision, dimension(:), pointer :: cvara_totwt, cpro_tempc, cvara_ntdrp

save ideb
data ideb/0/

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! --- Radiative fluxes are computed at the first call
!     and then for the frequency nfatr1

if (mod(ntcabs,nfatr1).eq.0.or.ideb.eq.0) then

  allocate(temray(kmx), qvray(kmx), qlray(kmx), ncray(kmx))
  allocate(fneray(kmx), romray(kmx), preray(kmx))
  allocate(zproj(kmx), ttvert(kmx*nvert), qvvert(kmx*nvert), romvert(kmx*nvert))
  allocate(qwvert(kmx*nvert), qlvert(kmx*nvert), ncvert(kmx*nvert))
  allocate(fnvert(kmx*nvert))
  allocate(aeroso(kmx))
  allocate(coords(3,kmx,nvert))
  allocate(cressm(kmx*nvert), interp(kmx*nvert))
  allocate(infrad(3*kmx*nvert))

  ideb = 1

  heuray = float(shour) + float(smin)/60.d0+ssec/3600.d0                       &
         + (ntcabs-1)*dtref/3600.d0

  if (ntcabs.le.2.or.isuite.eq.1) then
    ico2 = 1
  else
    ico2 = 0
  endif

  ! --- Initialization:
  do ii = 1, nvert
    soilvert(ii)%foir = 0.d0
    soilvert(ii)%fos  = 0.d0
  enddo

  do k = 2, kvert
    zray(k) = 0.d0
    preray(k) = 0.d0
    temray(k) = 0.d0
    qvray(k) = 0.d0
    romray(k) = 0.d0
    qlray (k) = 0.d0
    ncray (k) = 0.d0
    fneray(k) = 0.d0
    aeroso(k) = 0.d0
  enddo

  call field_get_val_s(icrom, crom)
  call field_get_val_s(itempc, cpro_tempc)

  if (ippmod(iatmos).eq.2) then
    call field_get_val_s(iliqwt, cpro_pcliq)

    call field_get_val_prev_s(ivarfl(isca(itotwt)), cvara_totwt)
    call field_get_val_prev_s(ivarfl(isca(intdrp)), cvara_ntdrp)
  endif

  !=============================================================================
  ! 2.  Computing long-wave and short-wave radiative fluxes
  !=============================================================================

  do ii = 1, nvert

    xvert = xyvert(ii,1)
    yvert = xyvert(ii,2)

    do jj = 1, kmx
      coords(1,jj,ii) = xvert
      coords(2,jj,ii) = yvert
      coords(3,jj,ii) = zvert(jj)
    enddo
  enddo

  if (ntcabs.eq.1) then
    call grimap(igrid, nvert*kmx, coords)
  endif

  call gripol(igrid, cpro_tempc, ttvert)
  call gripol(igrid, crom, romvert)

  if (ippmod(iatmos).eq.2) then
    ! liquid content interpolation
    call gripol(igrid, cpro_pcliq, qlvert)
    ! total water content interpolation
    call gripol(igrid, cvara_totwt, qwvert)

    ! deduce vapor content interpolation
    do ii = 1, nvert
      do k = 1, kvert
        jj = (ii-1)*kmx + k
        qvvert(jj) = qwvert(jj)-qlvert(jj)
      enddo
    enddo

    call gripol(igrid, cvara_ntdrp, ncvert)
    call gripol(igrid, nebdia, fnvert)
  endif

  ! --- Loop on the vertical array:

  do ii = 1, nvert
    xvert = xyvert(ii,1)
    yvert = xyvert(ii,2)

    ! Soil constants
    albedo = soilvert(ii)%albedo
    emis = soilvert(ii)%emissi
    fos = soilvert(ii)%fos
    imer1 = 0

    ! 2.1 Profiles used for the computation of the radiative fluxes
    !--------------------------------------------------------------

    ! Soil variables
    zray(1)   = zvert(1)
    temray(1) = soilvert(ii)%ttsoil
    qvray(1)  = soilvert(ii)%totwat
    romray(1) = soilvert(ii)%density
    preray(1) = soilvert(ii)%pressure
    qlray(1)  = 0.d0
    ncray(1)  = 0.d0
    fneray(1) = 0.d0
    aeroso(1) = 0.d0

    ! Interpolation of temperature, humidity, density on the vertical
    ! The ref pressure profile is the one computed from the meteo profile
    if (ippmod(iatmos).eq.2.and.moddis.eq.2) then
      qlray(1) = qlvert(1 + (ii-1)*kmx)
      ncray(1)  = ncvert(1 + (ii-1)*kmx)
      fneray(1) = fnvert(1 + (ii-1)*kmx)
    endif

    do k = 2, kvert
      zray(k) = zvert(k)

      temray(k) = ttvert(k + (ii-1)*kmx)
      qvray(k)  = qvvert(k + (ii-1)*kmx)
      romray(k) = romvert(k + (ii-1)*kmx)

      ! default values
      ncray(k) = 0.d0
      qlray(k) = 0.d0
      fneray(k) = 0.d0

      if (ippmod(iatmos).eq.2.and.moddis.eq.2) then
        ncray(k)  = ncvert(k + (ii-1)*kmx)
        qlray(k)  = qlvert(k + (ii-1)*kmx)
        fneray(k) = fnvert(k + (ii-1)*kmx)
      endif

      if (imeteo.eq.0) then
        call atmstd(zray(k), preray(k), dum, dum)
      else
        call intprf(nbmetd, nbmetm, ztmet, tmmet, phmet, zray(k), ttcabs, &
                    preray(k))
      endif
    enddo

    ! --- Filling the additional levels
    kmray = kmx

    do k = kvert+1, kmray
      zray(k) = zvert(k)
      qlray(k) = 0.d0
      qvray(k) = 0.d0
      fneray(k) = 0.d0
      aeroso(k) = 0.d0
      ncray(k) = 0.d0

      ! initialize with standard atmosphere
      call atmstd(zray(k), preray(k), temray(k), romray(k))
    enddo

    call cs_user_atmo_1d_rad_prf(preray, temray, romray, qvray, qlray, &
                                 ncray, aeroso)

    ! --- Smoothing the temperature and humidity profile in the damping zone

    if (imeteo.eq.1) then
      ktamp = 6
      do k = kvert - ktamp+1, kmray
        call intprf(nbmaxt, nbmetm, ztmet, tmmet,                              &
                    ttmet, zray(k), ttcabs, temray(k))
        call intprf(nbmaxt, nbmetm, ztmet, tmmet, qvmet,                       &
                    zray(k), ttcabs, qvray(k))
      enddo

      icompt = 0
      do k = kvert,2,-1
        icompt = icompt+1
        if (icompt.le.6) then
          zrac = 2.d0*(zray(k) - zray(nbmett-ktamp + 3))                       &
               /(zray(nbmett) - zray(nbmett - ktamp))
          fpond = (1.d0 + tanh(zrac))/2.d0
          temray(k) = ttvert(k + (ii-1)*kmx)*(1.d0 - fpond) + temray(k)*fpond
          qvray(k) = qvvert(k + (ii-1)*kmx)*(1.d0 - fpond) + qvray(k)*fpond
          qlray(k) = qlvert(k + (ii-1)*kmx)*(1.d0 - fpond) + qlray(k)*fpond
        endif
      enddo
    endif

    ! --- Clipping the humidity

    do k = 1, kmray
      qvray(k) = max(5.d-4,qvray(k))
    enddo

    ! --- Computing pressure and density according to theta and qv profiles

    do k = kvert-ktamp+1, kmray
      tmoy = 0.5d0*(temray(k-1)+temray(k)) + tkelvi
      rhum = rair*(1.d0+(rvsra-1.d0)*qvray(k))
      rap = -abs(gz)*(zray(k)-zray(k-1))/rhum/tmoy
      preray(k) = preray(k-1)*exp(rap)
      romray(k) = preray(k)/(temray(k) + tkelvi)/rhum
    enddo

    ! 2.2 Computing the radiative fluxes for the vertical
    !-----------------------------------------------------

    k1 = 1

    ! --- Long-wave: InfraRed
    call rayir(ii, k1, kmray, ico2, emis,                         &
               tauzq, tauz, tausup, zq,                           &
               acinfe, dacinfe, aco2, daco2, aco2s, daco2s,       &
               acsup, dacsup, acsups, dacsups,                    &
               zray, temray, qvray,                               &
               qlray, fneray, romray, preray, aeroso,             &
               foir, rayi(:,ii), ncray)

    ! --- Short-wave: Sun
    call rayso(ii, k1, kmray, heuray, imer1, albedo,             &
               tauzq, tauz, tausup, zq,                          &
               zray,                                             &
               qvray, qlray, fneray, romray, preray, aeroso,     &
               fos, rayst(:,ii), ncray)

    soilvert(ii)%fos = fos
    soilvert(ii)%foir = foir
  enddo

  do ii = 1, kmx*nvert
    cressm(ii) = 1
    interp(ii) = 1
    infrad(3*(ii-1) + 1) = 1.d0/8500.d0 ! horizontal(x) cressman radius
    infrad(3*(ii-1) + 2) = 1.d0/8500.d0 ! horizontal(y) cressman radius
    infrad(3*(ii-1) + 3) = 4.d0/1.d0    ! vertical(z) cressman radius
  enddo

  ! Map Infra Red 1D (rayi) on the structure idrayi
  call mesmap (idrayi, kmx*nvert, rayi, coords,            &
               cressm, interp, infrad)

  ! Map Sun 1D (rayst) on the structure idrayst
  call mesmap (idrayst, kmx*nvert, rayst, coords,          &
               cressm, interp, infrad)

  deallocate(temray, qvray, qlray)
  deallocate(fneray, romray, preray)
  deallocate(zproj, ttvert, qvvert, romvert)
  deallocate(aeroso)
  deallocate(cressm)
  deallocate(interp)
  deallocate(coords, infrad)

endif
end subroutine atr1vf
