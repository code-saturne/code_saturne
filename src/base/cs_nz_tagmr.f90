!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> cs_user_boundary_mass_source_terms.

module cs_nz_tagmr

  !=============================================================================

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
  integer, save :: znmurx

  !> \anchor znmur
  !> number of discretized points associated to the (ii)th face
  !> with the 1-D thermal model coupled with condensation
  integer, allocatable, dimension(:) :: znmur

  !> \anchor zdxp
  !> space step associated to the spatial discretization of the 1-D thermal model
  !> coupled with condensation model
  double precision, allocatable, dimension(:,:) :: zdxp

  !> \anchor ztheta-scheme of the 1-D thermal model
  !>    - 0 : explicit scheme
  !>    - 1 : implicit scheme
  double precision, allocatable, dimension(:) :: ztheta

  !> \anchor zdxmin
  !> the minimal space step of 1-D thermal model
  !> by default equal to 0 with a homogeneus space step.
  !> this numerical parameter is used to impose a geometric progression
  !> ratio of the mesh refinement associated to (ii)th face with
  !> the 1-D thermal model.
  double precision, allocatable, dimension(:) :: zdxmin

  !> \anchor zepais
  !> the wall thickness associated to the (ii)th face with 1-D thermal module
  double precision, allocatable, dimension(:) :: zepais

  !> \}

  !----------------------------------------------------------------------------
  ! Physical parameters to solve the 1-D thermal conduction equation
  !----------------------------------------------------------------------------

  !> \defgroup physical_properties of 1-D thermal conduction Equation

  !> \addtogroup physcial_properties of each zone
  !> \{

  !> \anchor zxrob
  !> the concrete density associated to solid material
  double precision, allocatable, dimension(:) :: zrob

  !> \anchor zcondb
  !> the concrete conductivity coefficient associated to solid material
  double precision, allocatable, dimension(:) :: zcondb

  !> \anchor zcpb
  !> the concrete specific heat coefficient associated to solid material
  double precision, allocatable, dimension(:) :: zcpb

  !> \anchor zhext
  !> the exterior exchange coefficient associated to solid material
  double precision, allocatable, dimension(:) :: zhext

  !> \anchor ztext
  !> the exterior temperature associated to solid material
  double precision, allocatable, dimension(:) :: ztext

  !> \anchor ztpar0
  !> the initial temperature associated to solid material
  double precision, allocatable, dimension(:) :: ztpar0

  !> \anchor ztmur
  !> the wall temperature computed with the 1-D thermal model
  !> associated to concrete solid material
  double precision, dimension(:,:),  allocatable :: ztmur

  !> \}

  !> \}

contains

  !=============================================================================

  subroutine init_nz_tagmr

    use cs_nz_condensation, only:nzones

    if (nzones.lt.1) nzones = 1

    allocate(znmur (nzones))
    allocate(ztheta(nzones))
    allocate(zdxmin(nzones))
    allocate(zepais(nzones))
    allocate(zrob  (nzones))
    allocate(zcondb(nzones))
    allocate(zcpb  (nzones))
    allocate(zhext (nzones))
    allocate(ztext (nzones))
    allocate(ztpar0(nzones))

    !---> Array initialization

    znmur(:)  = 0
    ztheta(:) = 0.d0
    zdxmin(:) = 0.d0
    zepais(:) = 0.d0
    zrob(:)   = 0.d0
    zcondb(:) = 0.d0
    zcpb(:)   = 0.d0
    zhext(:)  = 0.d0
    ztext(:)  = 0.d0
    ztpar0(:) = 0.d0

  end subroutine init_nz_tagmr

  !=============================================================================

  subroutine finalize_nz_tagmr

    deallocate(znmur )
    deallocate(ztheta)
    deallocate(zdxmin)
    deallocate(zepais)
    deallocate(zrob  )
    deallocate(zcondb)
    deallocate(zcpb  )
    deallocate(zhext )
    deallocate(ztext )
    deallocate(ztpar0)

  end subroutine finalize_nz_tagmr

  !=============================================================================

  subroutine init_nz_mesh_tagmr

    use optcal
    use pointe
    use cs_tagmr
    use cs_nz_condensation
    use parall

    implicit none

    ! Local variables

    integer  iiii, iz

    ! Copy single-zone to multi-zone formulation for compatibility if needed.

    if (znmur(1).eq.0) then
      nztag1d = itag1d
      do iiii = 1, nfbpcd
        iz = izzftcd(iiii)
        izcophc(iz) = icophc
        izcophg(iz) = icophg
        iztag1d(iz) = itag1d
        znmur (iz)  = nmur
        ztheta(iz)  = theta
        zdxmin(iz)  = dxmin
        zepais(iz)  = epais
        ztpar0(iz)  = tpar0
      enddo
    else
      ! Calcul du max des nmurs (pour les fichiers suite)
      nztag1d = 0
      do iz = 1, nzones
        do iiii = 1, nfbpcd
          if (izzftcd(iiii).eq.iz.and.iztag1d(iz).eq.1.and.iiii.gt.0) then
            nztag1d = max(iztag1d(iz), nztag1d)
          endif
        enddo
      enddo
      if (irangp.ge.0) call parcmx(nztag1d)
    endif

    ! the Condensation model coupled with a 1-D thermal model
    ! requires the 1-D mesh generation and temperature initialization
    if (nztag1d.eq.1) then

      if (nzones.eq.1) then
        znmurx = nmur
      else
        ! Calcul du max des nmurs (pour les fichiers suite)
        znmurx = 0
        do iz = 1, nzones
          znmurx = max(znmur(iz), znmurx)
        enddo
        if (irangp.ge.0) call parcmx(znmurx)
      endif

      ! Allocate zone data

      allocate(zdxp(nzones,znmurx))
      allocate(ztmur(nfbpcd,znmurx))

      zdxp(:,:) = 0.d0
      ztmur(:,:) = 0.d0

      !1-D mesh generated and temperature initialization
      call cs_mesh_tagmr(nfbpcd, izzftcd)

    endif

  end subroutine init_nz_mesh_tagmr

  !=============================================================================

  subroutine finalize_nz_mesh_tagmr

    deallocate(zdxp)
    deallocate(ztmur)

  end subroutine finalize_nz_mesh_tagmr

  !=============================================================================

end module cs_nz_tagmr
