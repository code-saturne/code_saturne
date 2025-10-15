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
!> \file atr1vf.f90
!>
!> \brief Compute radiative fluxes for the atmospheric model.
!> Computes the source term for scalar equations from radiative forcing
!> (UV and IR radiative fluxes) with a 1D scheme.
!-------------------------------------------------------------------------------

subroutine atr1vf () &
  bind(C, name="cs_f_atr1vf")

!===============================================================================
! Module files
!===============================================================================

use optcal
use cstphy
use cstnum, only:pi
use parall
use ppincl
use mesh
use field
use atincl
use cs_c_bindings
use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

procedure() :: grimap, gripol, rayir, rayso, mesmap

! Local variables
integer k, ii, jj
integer k1
integer ifac, isol
integer imer1
integer ideb, icompt
integer ktamp
integer nfmodsol, nbrsol

double precision heuray, albedo, emis, foir, fos
double precision xvert, yvert
double precision surf_zone
double precision zrac,fpond,rap,tmoy,rhum,dum
double precision soil_mean_albedo
double precision soil_mean_emissi
double precision soil_mean_density
double precision soil_mean_ttsoil
double precision soil_mean_totwat

integer, allocatable :: cressm(:), interp(:)
integer, dimension(:), pointer :: elt_ids

double precision, allocatable, dimension(:) :: temray, qvray, qlray, ncray
double precision, allocatable, dimension(:) :: fneray, romray, preray
double precision, allocatable, dimension(:) :: zproj, ttvert, romvert
double precision, allocatable, dimension(:) :: aeroso, infrad
double precision, allocatable, dimension(:,:,:) :: coords(:,:,:)
double precision, dimension(:), pointer :: crom, cpro_pcliq
double precision, dimension(:), pointer :: cvara_totwt, cpro_tempc, cvara_ntdrp
double precision, dimension(:), pointer :: nebdia
double precision, pointer, dimension(:)   :: bvar_tempp
double precision, pointer, dimension(:)   :: bvar_total_water
double precision, pointer, dimension(:)   :: bpro_albedo ! all boundary faces
double precision, pointer, dimension(:)   :: bpro_emissi ! all boundary faces
integer(c_int), dimension(2) :: dim_kmx_nvert
type(c_ptr) :: c_qwvert
type(c_ptr) :: c_qlvert
type(c_ptr) :: c_qvvert
type(c_ptr) :: c_ncvert
type(c_ptr) :: c_fnvert
type(c_ptr) :: c_aevert
double precision, dimension(:,:), pointer :: qwvert
double precision, dimension(:,:), pointer :: qlvert
double precision, dimension(:,:), pointer :: qvvert
double precision, dimension(:,:), pointer :: ncvert
double precision, dimension(:,:), pointer :: fnvert
double precision, dimension(:,:), pointer :: aevert

save ideb
data ideb/0/

interface

  subroutine cs_f_atmo_rad_1d_arrays_get_pointers(p_qwvert, p_qlvert, &
                                                  p_qvvert, p_ncvert, &
                                                  p_fnvert, p_aevert) &
    bind(C, name='cs_f_atmo_rad_1d_arrays_get_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: p_qwvert, p_qlvert
    type(c_ptr), intent(out) :: p_qvvert, p_ncvert
    type(c_ptr), intent(out) :: p_fnvert, p_aevert
  end subroutine cs_f_atmo_rad_1d_arrays_get_pointers

  subroutine cs_user_atmo_1d_rad_prf(preray, temray, romray, qvray, &
                                     qlray,  ncray, aeroso) &
    bind(C, name='cs_user_atmo_1d_rad_prf')
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=c_double), dimension(*), intent(inout) :: preray, temray
    real(kind=c_double), dimension(*), intent(inout) :: romray, qvray
    real(kind=c_double), dimension(*), intent(inout) :: qlray,  ncray, aeroso
  end subroutine cs_user_atmo_1d_rad_prf

end interface

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! --- Radiative fluxes are computed at the first call
!     and then for the frequency nfatr1

if (mod(ntcabs,nfatr1).eq.0.or.ideb.eq.0) then

  allocate(temray(kmx), qvray(kmx), qlray(kmx), ncray(kmx))
  allocate(fneray(kmx), romray(kmx), preray(kmx))
  allocate(zproj(kmx), ttvert(kmx*nvert), romvert(kmx*nvert))
  allocate(aeroso(kmx))

  ! Get pointers from C
  call cs_f_atmo_rad_1d_arrays_get_pointers(c_qwvert,  &
                                            c_qlvert,  &
                                            c_qvvert,  &
                                            c_ncvert,  &
                                            c_fnvert,  &
                                            c_aevert)

  dim_kmx_nvert(1) = kmx
  dim_kmx_nvert(2) = nvert
  call c_f_pointer(c_qwvert, qwvert, [dim_kmx_nvert])
  call c_f_pointer(c_qlvert, qlvert, [dim_kmx_nvert])
  call c_f_pointer(c_qvvert, qvvert, [dim_kmx_nvert])
  call c_f_pointer(c_ncvert, ncvert, [dim_kmx_nvert])
  call c_f_pointer(c_fnvert, fnvert, [dim_kmx_nvert])
  call c_f_pointer(c_aevert, aevert, [dim_kmx_nvert])

  allocate(coords(3,kmx+1,nvert))
  allocate(cressm(kmx*nvert), interp(kmx*nvert))
  allocate(infrad(3*kmx*nvert))

  ideb = 1

  heuray = float(shour) + float(smin)/60.d0+ssec/3600.d0

  if (idtvar.eq.0 .or. idtvar.eq.1) heuray = heuray + ttcabs/3600.d0

  ! --- Initialization:
  do k = 2, kvert
    preray(k) = 0.d0
    temray(k) = 0.d0
    qvray(k) = 0.d0
    romray(k) = 0.d0
    qlray (k) = 0.d0
    ncray (k) = 0.d0
    fneray(k) = 0.d0
  enddo

  call field_get_val_s_by_name("density", crom)
  call field_get_val_s_by_name("real_temperature", cpro_tempc)

  if (ippmod(iatmos).eq.2) then
    call field_get_val_s_by_name("liquid_water", cpro_pcliq)

    call field_get_val_prev_s_by_name("ym_water", cvara_totwt)
    call field_get_val_prev_s_by_name("number_of_droplets", cvara_ntdrp)
    call field_get_val_s_by_name("nebulosity_diag", nebdia)
  endif

  !=============================================================================
  ! 2.  Computing long-wave and short-wave radiative fluxes
  !=============================================================================

  ! Index of the bottom level
  k1 = 1

  do ii = 1, nvert ! (ixj) index

    xvert = xyvert(ii,1)
    yvert = xyvert(ii,2)

    ! Addition of one level for solar radiation
    ! TODO merge with 44km
    zq(kmx+1) = 16000.d0

    ! coords are levels (faces in 3D) whereas zray is slice (cells in 3D)
    do k = 1, kmx+1
      coords(1,k,ii) = xvert
      coords(2,k,ii) = yvert
      coords(3,k,ii) = zq(k)
      if (k.le.kmx) then
        zray(k) = 0.5d0 * (zq(k) + zq(k+1))
      endif
    enddo
  enddo

  if (ntcabs.eq.1) then
    call grimap(igrid, nvert*kmx, coords) ! face grid
  endif

  ! Grid interpolation is refurbished to get a P0 interpolation
  ! and not a point-point interpolation
  ! (i.e. we need to compute the mean of all code_saturne cells for a layer)
  ! A way could be to tag all nodes of the mesh is they belong to a ijk cell
  ! and then consider that a cell (or a face) belongs to a ijk if one node is
  ! in it (i.e. there is a non-zero intersection between the two).
  ! Then we can renormalise by the ratio "Volume(ijk)/SUM(cell_f_vol)"

  ! Interpolation from 3D to 1D in the fluid domain
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
        qvvert(k, ii) = qwvert(k, ii)-qlvert(k, ii)
      enddo
    enddo

    call gripol(igrid, cvara_ntdrp, ncvert)
    call gripol(igrid, nebdia, fnvert)
  endif

  ! Interpolate soil (P0 interpolation), same for all verticals for the moment

  if (iatsoil.ge.1) then

    call atmo_get_soil_zone(nfmodsol, nbrsol, elt_ids)

    ! Note: we use previous values of the soil to be coherent with
    ! previous datasetting and what have seen the fluid.
    call field_get_val_prev_s_by_name("soil_pot_temperature", bvar_tempp)
    call field_get_val_prev_s_by_name("soil_total_water", bvar_total_water)

    call field_get_val_s_by_name("boundary_albedo", bpro_albedo)
    call field_get_val_s_by_name("emissivity", bpro_emissi)

    soil_mean_albedo = 0.d0
    soil_mean_emissi = 0.d0

    soil_mean_density = 0.d0
    soil_mean_ttsoil  = 0.d0
    soil_mean_totwat  = 0.d0
    surf_zone = 0.d0

    do isol = 1, nfmodsol

      ifac = elt_ids(isol) + 1 ! C > Fortran

      ! Density: property of the boundary face:
      soil_mean_density = soil_mean_density + surfbn(ifac) * crom(ifabor(ifac))

      ! Potential temperature for consistency with the code before, TODO is it correct?
      soil_mean_ttsoil  = soil_mean_ttsoil  + surfbn(ifac) * (bvar_tempp(isol) - tkelvi)
      soil_mean_totwat  = soil_mean_totwat  + surfbn(ifac) * bvar_total_water(isol)

      soil_mean_albedo = soil_mean_albedo + surfbn(ifac) * bpro_albedo(ifac)
      soil_mean_emissi = soil_mean_emissi + surfbn(ifac) * bpro_emissi(ifac)

      ! Surface of the zone, could use it directly in C
      surf_zone = surf_zone + surfbn(ifac)
    enddo

    if (irangp.ge.0) then
      call parsom(soil_mean_density)
      call parsom(soil_mean_ttsoil )
      call parsom(soil_mean_totwat )
      call parsom(soil_mean_albedo )
      call parsom(soil_mean_emissi )
      call parsom(surf_zone)
    endif

    soil_mean_density = soil_mean_density / surf_zone
    soil_mean_ttsoil  = soil_mean_ttsoil  / surf_zone
    soil_mean_totwat  = soil_mean_totwat  / surf_zone
    soil_mean_albedo  = soil_mean_albedo  / surf_zone
    soil_mean_emissi  = soil_mean_emissi  / surf_zone
    ! For now, all verticals have the same value
    ! TODO: automatic treatment for pressure?
    do ii = 1, nvert
      soil_albedo (ii) = soil_mean_albedo
      soil_emissi (ii) = soil_mean_emissi
      soil_ttsoil (ii) = soil_mean_ttsoil
      soil_totwat (ii) = soil_mean_totwat
      soil_density(ii) = soil_mean_density
    enddo
  endif

  ! --- Loop on the vertical array:

  do ii = 1, nvert

    ! FIXME the x, y position plays no role...
    ! interpolation must be reviewed
    xvert = xyvert(ii,1)
    yvert = xyvert(ii,2)

    ! Soil constants
    albedo = soil_albedo(ii)
    emis   = soil_emissi(ii)

    imer1 = 0

    ! 2.1 Profiles used for the computation of the radiative fluxes
    !--------------------------------------------------------------

    ! Loop over the variable defined until the top of the full domain
    do k = 1, kmx

      aeroso(k) = aevert(k, ii)
      if (ippmod(iatmos).eq.2.and.moddis.eq.2) then
        qlray(k)  = qlvert(k, ii)
        ncray(k)  = ncvert(k, ii)
        fneray(k) = fnvert(k, ii)
      else
        ! default values
        ncray(k) = 0.d0
        qlray(k) = 0.d0
        fneray(k) = 0.d0
      endif

    enddo

    ! Soil variables
    !temray(1) = soil_ttsoil(ii)
    !qvray(1)  = soil_totwat(ii)
    !romray(1) = soil_density (ii)
    !preray(1) = soil_pressure(ii)

    do k = 1, kvert

      temray(k) = ttvert(k + (ii-1)*kmx)

      qvray(k)  = qvvert(k, ii)
      romray(k) = romvert(k + (ii-1)*kmx)

      if (imeteo.eq.0) then
        call atmstd(0.d0, &
                    p0,   &
                    t0,   &
                    zray(k), preray(k), dum, dum)
      else if (imeteo.eq.1) then
        call intprf(nbmetd, nbmetm, ztmet, tmmet, phmet, zray(k), ttcabs, &
                    preray(k))
      else
        !TODO would be more coherent with an averaging of "meteo_pressure"
        ! Field
        call atmstd(0.d0, &
                    p0,   &
                    t0,   &
                    zray(k), preray(k), dum, dum)
      endif
    enddo

    ! --- Filling the additional levels
    !TODO do it before
    do k = kvert+1, kmx

      ! Initialize with standard atmosphere
      ! above the domain
      call atmstd(zray(kvert), preray(kvert), temray(kvert)+tkelvi, &
                  zray(k), preray(k), temray(k), romray(k))
      ! Conversion Kelvin to Celsius
      temray(k) = temray(k) - tkelvi
    enddo

    call cs_user_atmo_1d_rad_prf(preray, temray, romray, qvray, qlray, &
                                 ncray, aeroso)

    ! --- Smoothing the temperature and humidity profile in the damping zone

    ktamp = 0
    if (imeteo.eq.1) then
      ktamp = min(6, kvert)
      do k = kvert - ktamp+1, kmx
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
          qvray(k) = qvvert(k, ii)*(1.d0 - fpond) + qvray(k)*fpond
          qlray(k) = qlvert(k, ii)*(1.d0 - fpond) + qlray(k)*fpond
        endif
      enddo
    endif

    ! --- Clipping the humidity

    do k = 1, kmx
      qvray(k) = max(0.d0,qvray(k))
      qlray(k) = max(0.d0,qlray(k))
    enddo

    ! --- Computing pressure and density according to temperature
    !     and qv profiles

    do k = kvert-ktamp+1, kmx
      tmoy = 0.5d0*(temray(k-1)+temray(k)) + tkelvi
      rhum = rair*(1.d0+(rvsra-1.d0)*qvray(k))
      rap = -abs(gz)*(zray(k)-zray(k-1))/rhum/tmoy
      preray(k) = preray(k-1)*exp(rap)
      romray(k) = preray(k)/(temray(k) + tkelvi)/rhum
    enddo

    ! 2.2 Computing the radiative fluxes for the vertical
    !-----------------------------------------------------

    ! --- Long-wave: InfraRed
    call rayir(ii, k1, kmx, emis,                                 &
               tauzq, tauz, tausup, zq,                           &
               zray, temray, qvray,                               &
               qlray, fneray, romray, preray, aeroso,             &
               soil_ttsoil(ii), soil_totwat(ii),                  &
               soil_density (ii), soil_pressure(ii),              &
               foir, rayi(:,ii), ncray)

    ! --- Short-wave: Sun
    call rayso(ii, k1, kmx, heuray, imer1, albedo,                &
               tauzq, tauz, tausup, zq,                           &
               zray,                                              &
               qvray, qlray, fneray, romray, preray, temray,      &
               fos, rayst(:,ii), ncray)

  enddo

  do ii = 1, kmx*nvert
    cressm(ii) = 1
    interp(ii) = 1
    infrad(3*(ii-1) + 1) = 1.d0/8500.d0 ! horizontal(x) Cressman radius
    infrad(3*(ii-1) + 2) = 1.d0/8500.d0 ! horizontal(y) Cressman radius
    infrad(3*(ii-1) + 3) = 4.d0/1.d0    ! vertical(z) Cressman radius
  enddo

  ! Map Infra Red 1D (rayi) on the structure idrayi
  call mesmap (idrayi, kmx*nvert, rayi, coords,            &
               cressm, interp, infrad)

  ! Map Sun 1D (rayst) on the structure idrayst
  call mesmap (idrayst, kmx*nvert, rayst, coords,          &
               cressm, interp, infrad)

  deallocate(temray, qvray, qlray)
  deallocate(fneray, romray, preray)
  deallocate(zproj, ttvert, romvert)
  deallocate(aeroso)
  deallocate(cressm)
  deallocate(interp)
  deallocate(coords, infrad)

endif
end subroutine atr1vf
