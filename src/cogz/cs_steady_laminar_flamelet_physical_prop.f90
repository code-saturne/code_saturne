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

!===============================================================================
! Function:
! ---------

!> \file cs_steady_laminar_flamelet_physical_prop.f90
!>
!> \brief Specific physic subroutine: diffusion flame.
!>
!> Calculation of mean thermophysic properties
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!______________________________________________________________________________!

subroutine cs_steady_laminar_flamelet_physical_prop()  &
  bind(C, name='cs_steady_laminar_flamelet_physical_prop')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppthch
use coincl
use ppincl
use radiat
use mesh
use field
use pointe, only:pmapper_double_r1
use parall
use period
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer           iel, ifcvsl, isc, iscal, iesp, ig, iprev
integer           viscls_counter, viscls_number

double precision  had, cmin, cmid, cmax

double precision, dimension(:), pointer :: cvar_fm, fp2m, cvar_progvar
double precision, dimension(:), pointer :: cpro_temp, cpro_progvar, cpro_omegac
double precision, dimension(:), pointer :: cvar_scalt, cpro_xr, cpro_totki
double precision, dimension(:), pointer :: cpro_rho, cpro_hrr
double precision, dimension(:), pointer :: cpro_viscl
double precision, dimension(:), pointer :: cpro_tem2

double precision, allocatable, dimension(:)  :: phim, rom_eos
double precision, allocatable, dimension(:,:)  :: rad_work

type(pmapper_double_r1), dimension(:), pointer :: cpro_kg, cpro_emi
type(pmapper_double_r1), dimension(:), pointer :: cpro_species
type(pmapper_double_r1), dimension(:), pointer :: cpro_viscls

character(len=64) :: f_name

logical          update_rad
integer, save :: ntcabs_prev = -1
logical, save :: is_vploop = .false.

interface

  subroutine max_mid_min_progvar(zm0, cmax, cmid, cmin)  &
    bind(C, name='cs_f_max_mid_min_progvar')
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(in)    :: zm0
    real(c_double), intent(inout) :: cmax, cmid, cmin
  end subroutine max_mid_min_progvar

  subroutine combustion_reconstruct_variance(iprev)  &
    bind(C, name='cs_f_combustion_reconstruct_variance')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int) :: iprev
  end subroutine combustion_reconstruct_variance

  subroutine cs_combustion_boundary_conditions_density()  &
    bind(C, name='cs_combustion_boundary_conditions_density')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_combustion_boundary_conditions_density

end interface

if (ntcabs_prev .eq. ntcabs) then
  is_vploop = .true.
else
  is_vploop = .false.
endif

ntcabs_prev = ntcabs

!===============================================================================
! 1. Map the pointers for fields
!===============================================================================

call field_get_val_s(ifm, cvar_fm)

if(mode_fp2m.eq.0) then
  call field_get_val_s(ifp2m, fp2m)
elseif(mode_fp2m.eq.1) then
  call field_get_val_s(irecvr, fp2m)
  iprev = 0
  call combustion_reconstruct_variance(iprev)
endif

if (ippmod(islfm).eq.1 .or. ippmod(islfm).eq.3)  then
  call field_get_val_s(ihm, cvar_scalt)
  call field_get_val_s(ixr, cpro_xr)
  do iel =1,ncel
    had = cvar_fm(iel) * hinfue + (1.d0-cvar_fm(iel))*hinoxy
    cpro_xr(iel) = max(-(cvar_scalt(iel) - Had),0.d0)
  enddo
else
  allocate(cpro_xr(ncelet))
  cpro_xr(:) = 0.d0
endif

if (ippmod(islfm).ge.2) then
  call field_get_val_s(ipvm, cvar_progvar)
  call field_get_val_s(iomgc, cpro_omegac)

  ! Clip the progress variable here, Not in cs_scalar_clipping
  do iel = 1, ncel
    call max_mid_min_progvar(cvar_fm(iel), cmax, cmid, cmin)
    cvar_progvar(iel) = min(cvar_progvar(iel), cmax)
  enddo
else
  call field_get_val_s(itotki, cpro_totki)
  call scalar_dissipation_rate
endif

call field_get_val_s(itemp , cpro_temp)

call field_get_val_s(iviscl, cpro_viscl)

call field_get_val_s(it2m,  cpro_tem2)

call field_get_val_s(iym(ngazgm), cpro_progvar)

call field_get_val_s(ihrr, cpro_hrr)

if (iirayo.eq.1 .and. mod(ntcabs, nt_rad_prp).eq.0) then
  update_rad = .true.
else
  update_rad = .false.
endif

!==============================================================================

if (is_vploop .eqv. .false.) then ! Outside the rho(Y)-v-p coupling

  !------------------------------------------------------------!
  ! Map pointers for multi-dimensional fields                  !
  ! to avoid repeatedly retrieving pointers in ncel loop below !
  !------------------------------------------------------------!

  ! Map species pointer
  allocate(cpro_species(ngazfl))
  do iesp = 1, ngazfl
    call field_get_val_s(iym(iesp), cpro_species(iesp)%p)
  enddo

  ! Map scalars' diffusivity pointer
  viscls_counter = 0
  do isc = 1, nscapp
    iscal = iscapp(isc)
    if (iscavr(iscal).le.0) then
      call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0) then
        viscls_counter = viscls_counter + 1
      endif
    endif
  enddo

  viscls_number = viscls_counter
  allocate(cpro_viscls(viscls_number))

  viscls_counter = 0
  do isc = 1, nscapp
    iscal = iscapp(isc)
    if (iscavr(iscal).le.0) then
      call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0) then
        viscls_counter = viscls_counter + 1
        call field_get_val_s(ifcvsl, cpro_viscls(viscls_counter)%p)
      endif
    endif
  enddo

  ! Map spectral absorption and emission pointers
  if (update_rad.eqv..TRUE.) then
    allocate(cpro_kg(nwsgg))
    allocate(cpro_emi(nwsgg))
    allocate(rad_work(2,nwsgg))

    do ig = 1, nwsgg
      write(f_name, '(a, i2.2)') 'spectral_absorption_coeff_', ig
      call field_get_val_s_by_name(f_name,  cpro_kg(ig)%p)

      write(f_name, '(a, i2.2)') 'spectral_emission_', ig
      call field_get_val_s_by_name(f_name, cpro_emi(ig)%p)
    enddo
  endif
  !------------------ End Mapping of pointers -----------------------

  allocate(phim(nlibvar))

  do iel = 1, ncel

    if (ippmod(islfm).lt.2) then
      call filtered_physical_prop(cvar_fm(iel),     fp2m(iel),                &
                                  cpro_totki(iel), cpro_xr(iel),              &
                                  phim, rad_work, update_rad)
    else
      call filtered_physical_prop_progvar(cvar_fm(iel), fp2m(iel),            &
                                          cvar_progvar(iel), cpro_xr(iel),    &
                                          phim, rad_work, update_rad)
    endif
    ! Affectation des variables dans les tableaux appropries

    ! Temperature : cpro_temp
    cpro_temp(iel) = phim(flamelet_temp)

    ! viscosite laminaire: cpro_viscl
    cpro_viscl(iel) = phim(flamelet_vis)

    ! Fraction of species
    if (ngazfl.ge.1) then
      do iesp = 1, ngazfl
        cpro_species(iesp)%p(iel) = phim(flamelet_species(iesp))
      enddo
    endif

    ! Diffusivite scalaire: cpro_viscls
    if (viscls_number .ge. 1) then
      do viscls_counter = 1, viscls_number
        cpro_viscls(viscls_counter)%p(iel) = phim(flamelet_dt)*phim(flamelet_rho)
      enddo
    endif

    if (flamelet_c.ne.-1) cpro_progvar(iel) = phim(flamelet_c)
    if (flamelet_hrr.ne.-1) cpro_hrr(iel) = phim(flamelet_hrr)
    if (flamelet_temp2.ne.-1) cpro_tem2(iel) = phim(flamelet_temp2)

    if (ippmod(islfm).ge.2) then
      cpro_omegac(iel) = phim(flamelet_omg_c)
    endif

    if (update_rad.eqv..TRUE.) then
      do ig = 1, nwsgg
        cpro_kg(ig)%p(iel) = rad_work(1,ig)
        cpro_emi(ig)%p(iel) = rad_work(2,ig)
      enddo
    endif
  enddo

  if(update_rad.eqv..TRUE.) then
    deallocate(rad_work)
    deallocate(cpro_kg, cpro_emi)
  endif

  deallocate(phim)
  deallocate(cpro_species)
  deallocate(cpro_viscls)

else  ! Inside the rho(Y)-v-p coupling

  call field_get_val_s(icrom, cpro_rho)

  allocate(rom_eos(ncelet))

  if (ippmod(islfm).lt.2) then
    do iel= 1, ncel
      call filtered_density(cvar_fm(iel), fp2m(iel), cpro_totki(iel), &
                            cpro_xr(iel), rom_eos(iel))
    enddo
  else
    do iel= 1, ncel
      call filtered_density_progvar(cvar_fm(iel), fp2m(iel),          &
                                    cvar_progvar(iel), cpro_xr(iel),  &
                                    rom_eos(iel))
    enddo
  endif

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rom_eos)
  endif

  if (itytur.eq.4) then
    call les_filter(1, rom_eos, cpro_rho)
  else
    do iel =1, ncel
      cpro_rho(iel) = rom_eos(iel)
    enddo
  endif

  deallocate(rom_eos)

  call cs_combustion_boundary_conditions_density

endif

if (ippmod(islfm).eq.0 .or. ippmod(islfm).eq.2) deallocate(cpro_xr)

return
end subroutine

!===============================================================================
! Function:
! ---------

!> \brief Specific physic subroutine: diffusion flame.
!>
!> Interpolation of mean thermophysic properties as a
!> function of fm, fp2m, ki, xr
!-------------------------------------------------------------------------------

subroutine filtered_physical_prop(zm0, zvar0, kim0, xrm0,         &
                                  phim, rad_work, update_rad)

!===============================================================================
! Module files
!===============================================================================

use coincl
use radiat

implicit none

integer :: dataIndex

double precision :: zm0, zvar0, kim0, xrm0, phim(nlibvar)
double precision :: weight_zm, weight_zvar, weight_kim, weight_xrm

double precision, dimension(:,:,:,:), allocatable :: phi_4
double precision, dimension(:,:,:), allocatable :: phi_3
double precision, dimension(:,:), allocatable :: phi_2
double precision, dimension(:), allocatable :: xData

! Work arrays for radiation
double precision, dimension(:,:,:,:,:), allocatable :: rad_work4
double precision, dimension(:,:,:,:), allocatable  :: rad_work3
double precision, dimension(:,:,:), allocatable    :: rad_work2
double precision  rad_work(2, nwsgg)

logical  update_rad

allocate(phi_4(nlibvar,nxr,nki,nzvar))
allocate(phi_3(nlibvar,nxr,nki))
allocate(phi_2(nlibvar,nxr))

if (update_rad.eqv..TRUE.) then
  allocate(rad_work4(2,nwsgg,nxr,nki,nzvar))
  allocate(rad_work3(2,nwsgg,nxr,nki))
  allocate(rad_work2(2,nwsgg,nxr))
endif

! Start with z_m
allocate(xData(nzm))
xData(:) = flamelet_library(flamelet_zm, 1, 1, 1, :)

if(xData(1).ge.zm0) then
  phi_4(:,:,:,:) = flamelet_library(:,:,:,:,1)
  if (update_rad.eqv..TRUE.) rad_work4(:,:,:,:,:) = radiation_library(:,:,:,:,:,1)
elseif(xData(size(xData)).le.zm0) then
  phi_4(:,:,:,:) = flamelet_library(:,:,:,:,nzm)
  if (update_rad.eqv..TRUE.) rad_work4(:,:,:,:,:) = radiation_library(:,:,:,:,:,nzm)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt.1 ) then
    weight_zm = (zm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_4(:,:,:,:) = (1.0-weight_zm)*flamelet_library(:,:,:,:,dataIndex-1) &
      + weight_zm*flamelet_library(:,:,:,:,dataIndex)

    ! Interpolation grandeurs radiatives
    if (update_rad.eqv..TRUE.) then
      rad_work4(:,:,:,:,:) = (1.0-weight_zm)*radiation_library(:,:,:,:,:,dataIndex-1) &
        + weight_zm*radiation_library(:,:,:,:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)

! Then interpolate over z_var
allocate(xData(nzvar))
xData(:) = phi_4(flamelet_zvar,1,1,:)

if(xData(1).ge.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,1)
  if (update_rad.eqv..TRUE.) rad_work3(:,:,:,:) = rad_work4(:,:,:,:,1)
elseif(xData(size(xData)).le.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,nzvar)
  if (update_rad.eqv..TRUE.) rad_work3(:,:,:,:) = rad_work4(:,:,:,:,nzvar)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zvar0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt. 1) then
    weight_zvar = (zvar0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_3(:,:,:) = (1.0-weight_zvar)*phi_4(:,:,:,dataIndex-1) &
      + weight_zvar*phi_4(:,:,:,dataIndex)
    ! Interpolation grandeurs radiatives
    if (update_rad.eqv..TRUE.) then
      rad_work3(:,:,:,:) = (1.0-weight_zvar)*rad_work4(:,:,:,:,dataIndex-1) &
        + weight_zvar*rad_work4(:,:,:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)

! Then interpolate over Ki_m
allocate(xData(nki))
xData(:) = phi_3(flamelet_ki,1,:)

if(xData(1).ge.kim0) then
  phi_2(:,:) = phi_3(:,:,1)
  if (update_rad.eqv..true.) rad_work2(:,:,:) = rad_work3(:,:,:,1)
elseif(xData(size(xData)).le.kim0) then
  phi_2(:,:) = phi_3(:,:,nki)
  if (update_rad.eqv..true.) rad_work2(:,:,:) = rad_work3(:,:,:,nki)
else

  dataIndex = 1
  do while(xData(dataIndex).lt. kim0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt. 1) then
    weight_kim = (kim0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_2(:,:) = (1.0-weight_kim)*phi_3(:,:,dataIndex-1) &
      + weight_kim*phi_3(:,:,dataIndex)
    ! Interpolation grandeurs radiatives
    if (update_rad.eqv..TRUE.) then
      rad_work2(:,:,:) = (1.0-weight_kim)*rad_work3(:,:,:,dataIndex-1) &
        + weight_kim*rad_work3(:,:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)

! Then interpolate over XR_m
allocate(xData(nxr))
xData(:) = phi_2(flamelet_xr,:)

if(xData(1).ge.xrm0) then
  phim(:) = phi_2(:,1)
  if (update_rad.eqv..TRUE.) rad_work(:,:) = rad_work2(:,:,1)
elseif(xData(size(xData)).le.xrm0) then
  phim(:) = phi_2(:,nxr)
  if (update_rad.eqv..TRUE.) rad_work(:,:) = rad_work2(:,:,nxr)
else

  dataIndex = 1
  do while(xData(dataIndex).lt. xrm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex.gt.1) then
    weight_xrm = (xrm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phim(:) = (1.0-weight_xrm)*phi_2(:,dataIndex-1) &
      + weight_xrm*phi_2(:,dataIndex)
    if (update_rad.eqv..TRUE.) then
      rad_work(:,:) = (1.0-weight_xrm)*rad_work2(:,:,dataIndex-1) &
        + weight_xrm*rad_work2(:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)
deallocate(phi_4,phi_3,phi_2)
if (update_rad.eqv..TRUE.) deallocate(rad_work4, rad_work3, rad_work2)

end subroutine

!===============================================================================
! Function:
! ---------

!> \brief Specific physic subroutine: diffusion flame.
!>
!> Interpolation of mean density as a function of
!> fm, fp2m, ki, xr
!-------------------------------------------------------------------------------

subroutine filtered_density(zm0, zvar0, kim0, xrm0, rhom)

!===============================================================================
! Module files
!===============================================================================

use coincl

implicit none

double precision :: zm0, zvar0, kim0, xrm0, rhom
double precision :: weight_zm, weight_zvar, weight_kim, weight_xrm

double precision, dimension(:,:,:,:), allocatable :: phi_4
double precision, dimension(:,:,:), allocatable :: phi_3
double precision, dimension(:,:), allocatable :: phi_2
double precision, dimension(:), allocatable :: xData

integer :: dataIndex, n_var_local = 5
integer :: var_sel(5)

var_sel = (/flamelet_zm, flamelet_zvar, flamelet_xr,          &
            flamelet_rho, flamelet_ki/)

allocate(phi_4(n_var_local,nxr,nki,nzvar))
allocate(phi_3(n_var_local,nxr,nki))
allocate(phi_2(n_var_local,nxr))

! Start with z_m
allocate(xData(nzm))
xData(:) = flamelet_library(flamelet_zm,1,1,1,:)

if(xData(1).ge.zm0) then
  phi_4(:,:,:,:) = flamelet_library(var_sel,:,:,:,1)
elseif(xData(size(xData)).le.zm0) then
  phi_4(:,:,:,:) = flamelet_library(var_sel,:,:,:,nzm)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt.1 ) then
    weight_zm = (zm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_4(:,:,:,:) = (1.0-weight_zm)*flamelet_library(var_sel,:,:,:,dataIndex-1) &
      + weight_zm*flamelet_library(var_sel,:,:,:,dataIndex)
  endif
endif

deallocate(xData)

! Then interpolate over z_var
allocate(xData(nzvar))
xData(:) = phi_4(2,1,1,:)

if(xData(1).ge.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,1)
elseif(xData(size(xData)).le.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,nzvar)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zvar0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt. 1) then
    weight_zvar = (zvar0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_3(:,:,:) = (1.0-weight_zvar)*phi_4(:,:,:,dataIndex-1) &
      + weight_zvar*phi_4(:,:,:,dataIndex)
  endif
endif

deallocate(xData)

! Then interpolate over Ki_m
allocate(xData(nki))
xData(:) = phi_3(5,1,:)

if(xData(1).ge.kim0) then
  phi_2(:,:) = phi_3(:,:,1)
elseif(xData(size(xData)).le.kim0) then
  phi_2(:,:) = phi_3(:,:,nki)
else

  dataIndex = 1
  do while(xData(dataIndex).lt. kim0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt. 1) then
    weight_kim = (kim0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_2(:,:) = (1.0-weight_kim)*phi_3(:,:,dataIndex-1) &
      + weight_kim*phi_3(:,:,dataIndex)
  endif
endif

deallocate(xData)

! Then interpolate over XR_m
allocate(xData(nxr))
xData(:) = phi_2(3,:)

if(xData(1).ge.xrm0) then
  rhom = phi_2(4,1)
elseif(xData(size(xData)).le.xrm0) then
  rhom = phi_2(4,nxr)
else

  dataIndex = 1
  do while(xData(dataIndex).lt. xrm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex.gt.1) then
    weight_xrm = (xrm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    rhom = (1.0-weight_xrm)*phi_2(4,dataIndex-1) &
      + weight_xrm*phi_2(4,dataIndex)
  endif
endif

deallocate(xData)
deallocate(phi_4,phi_3,phi_2)

end subroutine

!===============================================================================
! Function:
! ---------

!> \brief Specific physic subroutine: diffusion flame.
!>
!> Interpolation of mean thermophysic properties as a
!> function of fm, fp2m, c, xr
!-------------------------------------------------------------------------------

subroutine filtered_physical_prop_progvar(zm0, zvar0, cm0, xrm0,          &
                                          phim, rad_work, update_rad)

!===============================================================================
! Module files
!===============================================================================

use coincl
use radiat

implicit none

integer :: dataIndex, i_ki

double precision :: zm0, zvar0, cm0, xrm0, phim(nlibvar)
double precision :: weight_zm, weight_zvar, weight_cm, weight_xrm(nki)

double precision, dimension(:,:,:,:), allocatable :: phi_4
double precision, dimension(:,:,:), allocatable :: phi_3
double precision, dimension(:,:), allocatable :: phi_2, xData2D
double precision, dimension(:), allocatable :: xData

! Work arrays for radiation
double precision, dimension(:,:,:,:,:), allocatable :: rad_work4
double precision, dimension(:,:,:,:), allocatable  :: rad_work3
double precision, dimension(:,:,:), allocatable    :: rad_work2
double precision  rad_work(2, nwsgg)

logical update_rad

allocate(phi_4(nlibvar,nki,nxr,nzvar))
allocate(phi_3(nlibvar,nki,nxr))
allocate(phi_2(nlibvar,nki))

if (update_rad.eqv..TRUE.) then
  allocate(rad_work4(2,nwsgg,nki,nxr,nzvar))
  allocate(rad_work3(2,nwsgg,nki,nxr))
  allocate(rad_work2(2,nwsgg,nki))
endif

! Start with z_m
allocate(xData(nzm))
xData(:) = flamelet_library(flamelet_zm,1,1,1,:)

if(xData(1).ge.zm0) then
  phi_4(:,:,:,:) = flamelet_library(:,:,:,:,1)
  if (update_rad.eqv..TRUE.) rad_work4(:,:,:,:,:) = radiation_library(:,:,:,:,:,1)
elseif(xData(size(xData)).le.zm0) then
  phi_4(:,:,:,:) = flamelet_library(:,:,:,:,nzm)
  if (update_rad.eqv..TRUE.) rad_work4(:,:,:,:,:) = radiation_library(:,:,:,:,:,nzm)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo

  weight_zm = 0.d0
  if(dataIndex .gt.1 ) then
    weight_zm = (zm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_4(:,:,:,:) = (1.0d0-weight_zm)*flamelet_library(:,:,:,:,dataIndex-1) &
      + weight_zm*flamelet_library(:,:,:,:,dataIndex)

    ! Interpolation grandeurs radiatives
    if (update_rad.eqv..TRUE.) then
      rad_work4(:,:,:,:,:) = (1.0d0-weight_zm)*radiation_library(:,:,:,:,:,dataIndex-1) &
        + weight_zm*radiation_library(:,:,:,:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)

! Then interpolate over z_var
allocate(xData(nzvar))
xData(:) = phi_4(flamelet_zvar,1,1,:)

if(xData(1).ge.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,1)
  if (update_rad.eqv..TRUE.) rad_work3(:,:,:,:) = rad_work4(:,:,:,:,1)
elseif(xData(size(xData)).le.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,nzvar)
  if (update_rad.eqv..TRUE.) rad_work3(:,:,:,:) = rad_work4(:,:,:,:,nzvar)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zvar0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  weight_zvar = 0.d0
  if(dataIndex .gt. 1) then
    weight_zvar = (zvar0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_3(:,:,:) = (1.0d0-weight_zvar)*phi_4(:,:,:,dataIndex-1) &
      + weight_zvar*phi_4(:,:,:,dataIndex)

    ! Interpolation grandeurs radiatives
    if (update_rad.eqv..TRUE.) then
      rad_work3(:,:,:,:) = (1.0d0-weight_zvar)*rad_work4(:,:,:,:,dataIndex-1) &
        + weight_zvar*rad_work4(:,:,:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)

! Then interpolate over XR_m
allocate(xData2D(nki,nxr))
xData2D(:,:) = phi_3(4,:,:)

do i_ki = 1, nki

  if(xData2D(i_ki,1).ge.xrm0) then
    phi_2(:,i_ki) = phi_3(:,i_ki,1)
    if (update_rad.eqv..TRUE.) rad_work2(:,:,i_ki) = rad_work3(:,:,i_ki,1)
  elseif(xData2D(i_ki,nxr).le. xrm0) then
    phi_2(:,i_ki) = phi_3(:,i_ki,nxr)
    if (update_rad.eqv..TRUE.) rad_work2(:,:,i_ki) = rad_work3(:,:,i_ki,nxr)
  else

    dataIndex = 1
    do while(xData2D(i_ki,dataIndex) .lt. xrm0 &
        .and. dataIndex .lt. nxr)
      dataIndex = dataIndex + 1
    enddo

    if(dataIndex .gt. 1) then
      weight_xrm(i_ki) = (xrm0-xData2D(i_ki,dataIndex-1))/ &
        (xData2D(i_ki,dataIndex)-xData2D(i_ki,dataIndex-1))
      phi_2(:,i_ki) = (1.d0-weight_xrm(i_ki))*phi_3(:,i_ki,dataIndex-1) &
        + weight_xrm(i_ki)*phi_3(:,i_ki,dataIndex)
      if (update_rad.eqv..TRUE.) then
        rad_work2(:,:,i_ki) = (1.0d0-weight_xrm(i_ki))*rad_work3(:,:,i_ki,dataIndex-1) &
          + weight_xrm(i_ki)*rad_work3(:,:,i_ki,dataIndex)
      endif
    endif
  endif

enddo

deallocate(xData2D)

! Then interpolate over c_m
allocate(xData(nki))
xData(:) = phi_2(flamelet_c,:)

if(xData(1).ge.cm0) then
  phim(:) = phi_2(:,1)
  if (update_rad.eqv..true.) rad_work(:,:) = rad_work2(:,:,1)
elseif(xData(size(xData)).le.cm0) then
  phim(:) = phi_2(:,nki)
  phim(flamelet_omg_c) = 0.0d0
  if (update_rad.eqv..true.) rad_work(:,:) = rad_work2(:,:,nki)
else

  dataIndex = 1
  do while(xData(dataIndex).lt. cm0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  weight_cm = 0.d0
  if(dataIndex .gt. 1) then
    weight_cm = (cm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phim(:) = (1.0d0-weight_cm)*phi_2(:,dataIndex-1) &
      + weight_cm*phi_2(:,dataIndex)
    if (update_rad.eqv..TRUE.) then
      rad_work(:,:) = (1.0d0-weight_cm)*rad_work2(:,:,dataIndex-1) &
        + weight_cm*rad_work2(:,:,dataIndex)
    endif
  endif
endif

deallocate(xData)

deallocate(phi_4,phi_3,phi_2)
if (update_rad.eqv..TRUE.) deallocate(rad_work4, rad_work3, rad_work2)

end subroutine

!===============================================================================
! Function:
! ---------

!> \brief Specific physic subroutine: diffusion flame.
!>
!> Interpolation of mean density as a function of
!> fm, fp2m, c, xr
!-------------------------------------------------------------------------------

subroutine filtered_density_progvar(zm0, zvar0, cm0, xrm0, rhom)

!===============================================================================
! Module files
!===============================================================================

use coincl

implicit none

double precision, dimension(:), allocatable :: xData
double precision :: zm0, zvar0, cm0, xrm0, rhom

double precision, dimension(:,:,:,:), allocatable :: phi_4
double precision, dimension(:,:,:), allocatable :: phi_3
double precision, dimension(:,:), allocatable :: phi_2, xData2D
double precision :: weight_zm, weight_zvar, weight_cm, weight_xrm(nki)

integer :: dataIndex, i_ki, n_var_local = 5
integer :: var_sel(5)

var_sel = (/flamelet_zm, flamelet_zvar, flamelet_xr,          &
            flamelet_rho, flamelet_ki/)

allocate(phi_4(n_var_local,nki,nxr,nzvar))
allocate(phi_3(n_var_local,nki,nxr))
allocate(phi_2(n_var_local,nki))

! Start with z_m
allocate(xData(nzm))
xData(:) = flamelet_library(flamelet_zm,1,1,1,:)

if(xData(1).ge.zm0) then
  phi_4(:,:,:,:) = flamelet_library(var_sel,:,:,:,1)
elseif(xData(size(xData)).le.zm0) then
  phi_4(:,:,:,:) = flamelet_library(var_sel,:,:,:,nzm)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo

  weight_zm = 0.d0
  if(dataIndex .gt.1 ) then
    weight_zm = (zm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_4(:,:,:,:) = (1.0d0-weight_zm)*flamelet_library(var_sel,:,:,:,dataIndex-1) &
      + weight_zm*flamelet_library(var_sel,:,:,:,dataIndex)
  endif
endif

deallocate(xData)

! Then interpolate over z_var
allocate(xData(nzvar))
xData(:) = phi_4(2,1,1,:)

if(xData(1).ge.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,1)
elseif(xData(size(xData)).le.zvar0) then
  phi_3(:,:,:) = phi_4(:,:,:,nzvar)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zvar0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  weight_zvar = 0.d0
  if(dataIndex .gt. 1) then
    weight_zvar = (zvar0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    phi_3(:,:,:) = (1.0d0-weight_zvar)*phi_4(:,:,:,dataIndex-1) &
      + weight_zvar*phi_4(:,:,:,dataIndex)
  endif
endif

deallocate(xData)

! Then interpolate over XR_m
allocate(xData2D(nki,nxr))
xData2D(:,:) = phi_3(3,:,:)
do i_ki = 1, nki

  if(xData2D(i_ki,1).ge.xrm0) then
    phi_2(:,i_ki) = phi_3(:,i_ki,1)
  elseif(xData2D(i_ki,nxr).le. xrm0) then
    phi_2(:,i_ki) = phi_3(:,i_ki,nxr)
  else

    dataIndex = 1
    do while(xData2D(i_ki,dataIndex) .lt. xrm0 &
        .and. dataIndex .lt. nxr)
      dataIndex = dataIndex + 1
    enddo

    if(dataIndex .gt. 1) then
      weight_xrm(i_ki) = (xrm0-xData2D(i_ki,dataIndex-1))/ &
        (xData2D(i_ki,dataIndex)-xData2D(i_ki,dataIndex-1))
      phi_2(:,i_ki) = (1.d0-weight_xrm(i_ki))*phi_3(:,i_ki,dataIndex-1) &
        + weight_xrm(i_ki)*phi_3(:,i_ki,dataIndex)
    endif
  endif

enddo
deallocate(xData2D)

! Then interpolate over c_m
allocate(xData(nki))
xData(:) = phi_2(5,:)

if(xData(1).ge.cm0) then
  rhom = phi_2(4,1)
elseif(xData(size(xData)).le.cm0) then
  rhom = phi_2(4,nki)
else

  dataIndex = 1
  do while(xData(dataIndex).lt. cm0 &
      .and. dataIndex .lt. size(XData))
    dataIndex = dataIndex + 1
  enddo

  if(dataIndex .gt. 1) then
    weight_cm = (cm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    rhom = (1.0d0-weight_cm)*phi_2(4,dataIndex-1) &
      + weight_cm*phi_2(4,dataIndex)
  endif
endif

deallocate(xData)

deallocate(phi_4,phi_3,phi_2)

end subroutine

!===============================================================================
! Function:
! ---------

!> \brief Specific physic subroutine: diffusion flame.
!>
!> Calculation of total scalar dissipation rate
!-------------------------------------------------------------------------------

subroutine scalar_dissipation_rate ()

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppthch
use coincl
use ppincl
use field
use field_operator
use mesh
use parall
use period
use cs_c_bindings
use pointe

!====================================================================

implicit none

integer          iel, ifcvsl
integer          key_turb_diff, t_dif_id

double precision  dnum, delta_les

double precision, dimension(:), pointer :: cvar_fm, cpro_turb_diff, cpro_viscls

double precision, dimension(:), pointer :: fp2m, cpro_totki
double precision, dimension(:), pointer :: cpro_rhoa

!===================================================================

call field_get_val_prev_s(icrom, cpro_rhoa)
call field_get_val_s(ifm, cvar_fm)

if (mode_fp2m.eq.0) then
  call field_get_val_s(ifp2m, fp2m)
else if (mode_fp2m.eq.1) then
  call field_get_val_s(irecvr, fp2m)
endif

call field_get_key_int(ifm, kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

call field_get_val_s(itotki, cpro_totki)

call field_get_key_id("turbulent_diffusivity_id", key_turb_diff)
call field_get_key_int(ifm, key_turb_diff, t_dif_id)

if (t_dif_id .ge. 0) call field_get_val_s(t_dif_id, cpro_turb_diff)

if (iturb.eq.41) then
  do iel = 1, ncel
    delta_les = xlesfl *(ales*volume(iel))**bles

    dnum = coef_k*delta_les**2.d0*cpro_rhoa(iel)
    cpro_totki(iel) = (cpro_viscls(iel)+cpro_turb_diff(iel))/dnum * fp2m(iel)
  enddo
endif

return

end subroutine

!===============================================================================
! Function:
! ---------

!> \brief Specific physic subroutine: diffusion flame.
!>
!> Calculation of reconstructed variance of mixture fraction
!-------------------------------------------------------------------------------

subroutine combustion_reconstruct_variance(iprev) &
  bind(C, name='cs_f_combustion_reconstruct_variance')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppthch
use coincl
use ppincl
use field
use field_operator
use mesh
use parall
use period
use cs_c_bindings
use pointe

implicit none
integer          iel
integer(c_int) :: iprev

double precision, pointer, dimension(:) :: fm, fp2m, fsqm

call field_get_val_s(irecvr, fp2m)

if(iprev .eq. 0) then
  call field_get_val_s(ifm, fm)
  call field_get_val_s(ifsqm, fsqm)
else
  call field_get_val_prev_s(ifm, fm)
  call field_get_val_prev_s(ifsqm, fsqm)
endif

do iel = 1, ncel
  fp2m(iel) = fsqm(iel) - fm(iel)**2.d0
  fp2m(iel) = max(min(fp2m(iel), fm(iel)*(1.d0-fm(iel))),0.d0)
enddo

end subroutine


subroutine max_mid_min_progvar(zm0, cmax, cmid, cmin) &
    bind(C, name='cs_f_max_mid_min_progvar')

!===============================================================================
! Module files
!===============================================================================

use coincl

implicit none

real(c_double), intent(in)    :: zm0
real(c_double), intent(inout) :: cmax, cmid, cmin

double precision, dimension(:), allocatable :: xData
double precision :: weight_zm
integer :: dataIndex

! Start with z_m
allocate(xData(nzm))
xData(:) = flamelet_library(1,1,1,1,:)

if(xData(1).ge.zm0) then
  cmax = flamelet_library(flamelet_c,nki,1,1,1)
  cmid = flamelet_library(flamelet_c,ikimid,1,1,1)
  cmin = flamelet_library(flamelet_c,1,1,1,1)
elseif(xData(size(xData)).le.zm0) then
  cmax = flamelet_library(flamelet_c,nki,1,1,nzm)
  cmid = flamelet_library(flamelet_c,ikimid,1,1,nzm)
  cmin = flamelet_library(flamelet_c,1,1,1,nzm)
else

  dataIndex = 1
  do while(xData(dataIndex).lt.zm0 &
      .and. dataIndex .lt. size(xData))
    dataIndex = dataIndex + 1
  enddo
  if(dataIndex .gt.1 ) then
    weight_zm = (zm0 - xData(dataIndex-1))/(xData(dataIndex)-xData(dataIndex-1))
    cmax = (1.0d0-weight_zm)*flamelet_library(flamelet_c,nki,1,1,dataIndex-1) &
      + weight_zm*flamelet_library(flamelet_c,nki,1,1,dataIndex)

    cmid = (1.0d0-weight_zm)*flamelet_library(flamelet_c,ikimid,1,1,dataIndex-1) &
      + weight_zm*flamelet_library(flamelet_c,ikimid,1,1,dataIndex)

    cmin = (1.0d0-weight_zm)*flamelet_library(flamelet_c,1,1,1,dataIndex-1) &
      + weight_zm*flamelet_library(flamelet_c,1,1,1,dataIndex)

  endif
endif

deallocate(xData)

end subroutine max_mid_min_progvar



!end module steady_laminar
