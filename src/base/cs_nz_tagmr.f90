!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file cs_nz_tagmr.f90
!> Module for parameters options, numerical and physical properties of the
!> thermal 1D model for each specific zone with condensation on the wall.
!> The zones number is defined by the user with the subroutine :
!> cs_user_wall_condensation.

module cs_nz_tagmr

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================
  !> \defgroup cs_nz_tagmr Module for calculation options

  !> \addtogroup cs_nz_tagmr
  !> \{

  !----------------------------------------------------------------------------
  ! Numerical parameters to solve the 1-D thermal conduction equation
  !----------------------------------------------------------------------------

  !> \defgroup numerical_1d_nz parameters of the 1-D conduction equation

  !> \addtogroup numerical_1d_nz
  !> \{

  !> \anchor znmurx
  !> Maximal number of discretized points associated to the (ii)th face
  !> with the 1-D thermal model coupled with condensation
  integer(c_int), pointer, save :: znmurx

  !> \anchor znmur
  !> number of discretized points associated to the (ii)th face
  !> with the 1-D thermal model coupled with condensation
  integer, dimension(:), pointer, save :: znmur

  !> \anchor zdxp
  !> space step associated to the spatial discretization of the 1-D thermal model
  !> coupled with condensation model
  double precision, dimension(:,:), pointer, save :: zdxp

  !> \anchor ztheta-scheme of the 1-D thermal model
  !>    - 0 : explicit scheme
  !>    - 1 : implicit scheme
  double precision, dimension(:), pointer, save :: ztheta

  !> \anchor zdxmin
  !> the minimal space step of 1-D thermal model
  !> by default equal to 0 with a homogeneus space step.
  !> this numerical parameter is used to impose a geometric progression
  !> ratio of the mesh refinement associated to (ii)th face with
  !> the 1-D thermal model.
  double precision, dimension(:), pointer, save :: zdxmin

  !> \anchor zepais
  !> the wall thickness associated to the (ii)th face with 1-D thermal module
  double precision, dimension(:), pointer, save :: zepais

  !> \}

  !----------------------------------------------------------------------------
  ! Physical parameters to solve the 1-D thermal conduction equation
  !----------------------------------------------------------------------------

  !> \defgroup physical_properties of 1-D thermal conduction Equation

  !> \addtogroup physcial_properties of each zone
  !> \{

  !> \anchor zxrob
  !> the concrete density associated to solid material
  double precision, dimension(:), pointer, save :: zrob

  !> \anchor zcondb
  !> the concrete conductivity coefficient associated to solid material
  double precision, dimension(:), pointer, save :: zcondb

  !> \anchor zcpb
  !> the concrete specific heat coefficient associated to solid material
  double precision, dimension(:), pointer, save :: zcpb

  !> \anchor zhext
  !> the exterior exchange coefficient associated to solid material
  double precision, dimension(:), pointer, save :: zhext

  !> \anchor ztext
  !> the exterior temperature associated to solid material
  double precision, dimension(:), pointer, save :: ztext

  !> \anchor ztpar0
  !> the initial temperature associated to solid material
  double precision, dimension(:), pointer, save :: ztpar0

  !> \anchor ztmur
  !> the wall temperature computed with the 1-D thermal model
  !> associated to concrete solid material
  double precision, dimension(:,:), pointer, save :: ztmur

  !> \}

  !----------------------------------------------------------------------------
  ! Physical parameters to solve the 0-D thermal model of the metal mass
  !----------------------------------------------------------------------------

  !> \defgroup cs_tagms_physical_properties Physical properties of 0-D thermal model

  !> \addtogroup cs_tagms_physcial_properties
  !> \{

  !> \anchor  xem
  !> the wall thickness of the metal mass structures used by the 0-D thermal model
  double precision, dimension(:), pointer, save :: xem

  !> \anchor xro
  !> the density (kg.m-3) of the metal mass structures
  double precision, dimension(:), pointer, save :: xro_m

  !> \anchor xcond
  !> the conductivity coefficient (W.m-1.C-1) of the metal mass structures
  double precision, dimension(:), pointer, save :: xcond_m

  !> \anchor xcp
  !> the specific heat coefficient (J.kg-1.C-1) of the metal mass structures
  double precision, dimension(:), pointer, save :: xcp_m

  !> \anchor m_metal
  !> the metal mass (kg) of the metal mass structures
  double precision, dimension(:), pointer, save :: m_metal

  !> \anchor s_metal
  !> the exchange surface (m2) of the metal mass structures
  double precision, dimension(:), pointer, save :: s_metal

  !> \anchor v_metal
  !> the volume mesaure (m3) of the metal mass structures
  double precision, dimension(:), pointer, save :: v_metal

  !> \anchor tmet0
  !> the initial temperature of the metal mass structures
  double precision, dimension(:), pointer, save :: tmet0

  !> \anchor t_metal
  !> the wall temperature computed with the 0-D thermal model
  !> associated to the metal mass structure material
  double precision, dimension(:,:), pointer, save :: t_metal

  !> \}

  !> \}

interface

  subroutine cs_f_wall_condensation_1d_thermal_get_pointers(znmur, ztheta, zdxmin,&
                                                            zepais, zrob, zcondb,&
                                                            zcpb, zhext, ztext,&
                                                            ztpar0) &
    bind(C, name='cs_f_wall_condensation_1d_thermal_get_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: znmur, ztheta, zdxmin, zepais, zrob, zcondb
    type(c_ptr), intent(out) :: zcpb, zhext, ztext, ztpar0
  end subroutine cs_f_wall_condensation_1d_thermal_get_pointers

  subroutine cs_f_wall_condensation_1d_thermal_get_mesh_pointers(znmurx, zdxp, ztmur)&
    bind(C, name='cs_f_wall_condensation_1d_thermal_get_mesh_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: znmurx, zdxp, ztmur
  end subroutine cs_f_wall_condensation_1d_thermal_get_mesh_pointers

  subroutine cs_f_wall_condensation_1d_thermal_create(nzones) &
    bind(C, name='cs_wall_condensation_1d_thermal_create')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: nzones
  end subroutine cs_f_wall_condensation_1d_thermal_create

  subroutine cs_f_wall_condensation_0d_thermal_create(nvolumes, ncmast) &
    bind(C, name='cs_wall_condensation_0d_thermal_create')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: nvolumes, ncmast
  end subroutine cs_f_wall_condensation_0d_thermal_create

  subroutine cs_f_wall_condensation_1d_thermal_mesh_create(znmurx, nfbpcd, nzones) &
    bind(C, name='cs_wall_condensation_1d_thermal_mesh_create')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: znmurx, nfbpcd, nzones
  end subroutine cs_f_wall_condensation_1d_thermal_mesh_create

  subroutine cs_f_wall_condensation_0d_thermal_get_pointers(xem, t_metal,     &
                                                            xro_m, xcp_m,     &
                                                            xcond_m, m_metal, &
                                                            s_metal, v_metal, &
                                                            tmet0) &
    bind(C, name='cs_f_wall_condensation_0d_thermal_get_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: xem, t_metal, xro_m, xcp_m, xcond_m
    type(c_ptr), intent(out) :: m_metal, s_metal, v_metal, tmet0
  end subroutine cs_f_wall_condensation_0d_thermal_get_pointers

  subroutine cs_f_wall_condensation_0d_thermal_solve() &
    bind(C, name='cs_wall_condensation_0d_thermal_solve')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_f_wall_condensation_0d_thermal_solve

end interface


contains

  !=============================================================================

  subroutine init_nz_tagmr

    use, intrinsic :: iso_c_binding
    use cs_nz_condensation

    implicit none
    type(c_ptr) :: c_znmur, c_ztheta, c_zdxmin
    type(c_ptr) :: c_zepais, c_zrob, c_zcondb
    type(c_ptr) :: c_zcpb, c_zhext, c_ztext, c_ztpar0
    type(c_ptr) :: c_xem, c_t_metal, c_xro_m, c_xcp_m, c_xcond_m
    type(c_ptr) :: c_m_metal, c_s_metal, c_v_metal, c_tmet0

    if (nzones.lt.1) nzones = 1
    if (nvolumes.lt.1) nvolumes = 1

    call cs_f_wall_condensation_1d_thermal_create(nzones)
    call cs_f_wall_condensation_1d_thermal_get_pointers(c_znmur, c_ztheta, c_zdxmin,&
                                                        c_zepais, c_zrob, c_zcondb,&
                                                        c_zcpb, c_zhext, c_ztext,&
                                                        c_ztpar0)

    call c_f_pointer(c_znmur , znmur , [nzones])
    call c_f_pointer(c_ztheta, ztheta, [nzones])
    call c_f_pointer(c_zdxmin, zdxmin, [nzones])
    call c_f_pointer(c_zepais, zepais, [nzones])
    call c_f_pointer(c_zrob  , zrob  , [nzones])
    call c_f_pointer(c_zcondb, zcondb, [nzones])
    call c_f_pointer(c_zcpb  , zcpb  , [nzones])
    call c_f_pointer(c_zhext , zhext , [nzones])
    call c_f_pointer(c_ztext , ztext , [nzones])
    call c_f_pointer(c_ztpar0, ztpar0, [nzones])

    call cs_f_wall_condensation_0d_thermal_create(nvolumes, ncmast)

    call cs_f_wall_condensation_0d_thermal_get_pointers(c_xem, c_t_metal, c_xro_m, &
                                                        c_xcp_m, c_xcond_m,  &
                                                        c_m_metal, c_s_metal, &
                                                        c_v_metal, c_tmet0)

    call c_f_pointer(c_xem, xem, [nvolumes])
    call c_f_pointer(c_t_metal, t_metal, [2,nvolumes])
    call c_f_pointer(c_xro_m, xro_m, [nvolumes])
    call c_f_pointer(c_xcp_m, xcp_m, [nvolumes])
    call c_f_pointer(c_xcond_m, xcond_m, [nvolumes])
    call c_f_pointer(c_m_metal, m_metal, [nvolumes])
    call c_f_pointer(c_s_metal, s_metal, [nvolumes])
    call c_f_pointer(c_v_metal, v_metal, [nvolumes])
    call c_f_pointer(c_tmet0, tmet0, [nvolumes])

  end subroutine init_nz_tagmr

  !=============================================================================

  subroutine init_nz_mesh_tagmr

    use, intrinsic :: iso_c_binding
    use optcal
    use pointe
    use cs_nz_condensation
    use parall

    implicit none

    procedure() :: cs_mesh_tagmr

    ! Local variables

    integer  iiii, iz
    type(c_ptr) :: c_znmurx, c_zdxp, c_ztmur, c_dummy1, c_dummy2

    ! Copy single-zone to multi-zone formulation for compatibility if needed.
    call cs_f_wall_condensation_1d_thermal_get_mesh_pointers(c_znmurx, c_dummy1, &
                                                             c_dummy2)
    call c_f_pointer(c_znmurx, znmurx)

    ! Calcul du max des nmurs (pour les fichiers suite)
    nztag1d = 0
    do iz = 1, nzones
      do iiii = 1, nfbpcd
        if ((izzftcd(iiii)+1).eq.iz.and.iztag1d(iz).eq.1.and.iiii.gt.0) then
          nztag1d = max(iztag1d(iz), nztag1d)
        endif
      enddo
    enddo
    if (irangp.ge.0) call parcmx(nztag1d)

    ! the Condensation model coupled with a 1-D thermal model
    ! requires the 1-D mesh generation and temperature initialization
    if (nztag1d.eq.1) then

      ! Calcul du max des nmurs (pour les fichiers suite)
      znmurx = 0
      do iz = 1, nzones
        znmurx = max(znmur(iz), znmurx)
      enddo
      if (irangp.ge.0) call parcmx(znmurx)

      call cs_f_wall_condensation_1d_thermal_mesh_create(znmurx, nfbpcd, nzones)
      call cs_f_wall_condensation_1d_thermal_get_mesh_pointers(c_dummy1, c_zdxp, c_ztmur)
      call c_f_pointer(c_zdxp, zdxp, [nzones, znmurx])
      call c_f_pointer(c_ztmur, ztmur, [nfbpcd, znmurx])

      !1-D mesh generated and temperature initialization
      call cs_mesh_tagmr(nfbpcd, izzftcd)

    endif

  end subroutine init_nz_mesh_tagmr

  !=============================================================================

  subroutine finalize_nz_mesh_tagmr
    return
  end subroutine finalize_nz_mesh_tagmr

  !=============================================================================

end module cs_nz_tagmr
