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

!> \file cs_nz_condensation.f90
!> Module for parameters options and physical properties of the condensation
!> model using specific zones with different wall temperatures and material
!> properties.

module cs_nz_condensation

  !=============================================================================

  implicit none

  !=============================================================================
  !> \defgroup cs_nz_condensation Module for calculation options

  !> \addtogroup cs_nz_condensation
  !> \{

  !----------------------------------------------------------------------------
  ! Parameters to compute the condensation source terms
  !----------------------------------------------------------------------------

  !> \defgroup Parameters to choice the exchange coefficients type and the
  !> wall temperature used by the condensation model
  !> (with or without 1D thermal model to the wall).

  !> \addtogroup condensation and thermal exchange coefficients
  !> \{

  !> \anchor nzones
  !> number of the zones with a specific condensation source terms
  !> depending on the wall temperature and material properties.
  !> by default (nzones = 1) if the user does not specified different zones
  !> in the user subroutine \ref cs_user_boundary_mass_source_terms.f90 .
  integer, save :: nzones

  !> \anchor izzftcd
  !> list of the zones associated to the faces where a condensation source term
  !> is imposed. On each zone a specific wall temperature and material properties
  !> can be specified by the user.
  integer, allocatable, dimension(:) :: izzftcd

  !> \anchor izcophc
  !> choice the way to compute the exchange coefficient of the
  !> condensation source term used by the copain model.
  !>    - 1: the turbulent exchange coefficient of the flow
  !>    - 2: the exchange coefficient of the copain correlation
  !>    - 3: the maximal value between the two previous exchange coefficients
  integer, allocatable, dimension(:) :: izcophc

  !> \anchor izcophg
  !> choice the way to compute the thermal exchange coefficient associated
  !> to the heat transfer to the wall due to the condensation phenomenon.
  !>    - 2: the thermal exchange coefficient of the copain correlation
  !>    - 3: the maximal value between the current and previous thermal
  !>         exchange coefficient evaluated by the copain correlation
  integer, allocatable, dimension(:) :: izcophg

  !> \}

  !> \addtogroup Wall temperature model
  !> \{

  !> \anchor iztag1d
  !> choice the way to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the wall
  !>    - 1: the wall temperature is computed with a 1-D thermal model
  !>         with implicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         exchange coefficient evaluated by the copain correlation
  integer, allocatable, dimension(:) :: iztag1d

  !> \anchor nztag1d
  !> Indicate if the thermal 1D model of severe accident is usde
  !> to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the wall
  !>    - 1: the wall temperature is computed with a 1-D thermal model
  !>         with implicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         exchange coefficient evaluated by the copain correlation
  integer, save :: nztag1d

  !> Constant value of the wall temperature given by the user when
  !> the thermal 1D model is not activated for the condensation model with
  !> different zones specified in the user subroutine
  !> \ref cs_user_boundary_mass_source_terms.
  double precision, allocatable, dimension(:) :: ztpar

  !> \}
  !> \}

  !> \}

contains

  !=============================================================================

  subroutine init_nz_pcond

    use pointe, only:nfbpcd

    implicit none

    ! Local variables

    integer  ii

    allocate(izzftcd(nfbpcd))

    !---> Array initialization
    if (nzones.lt.1) then
      nzones = 1
      do ii = 1, nfbpcd
        izzftcd(ii) = 1
      enddo
    else
      izzftcd(:) = 0
    endif

    allocate(izcophc(nzones))
    allocate(izcophg(nzones))
    allocate(iztag1d(nzones))
    allocate(ztpar(nzones))

    izcophc(:) = 0
    izcophg(:) = 0
    iztag1d(:) = 0

    ztpar(:)  = -1.d0

  end subroutine init_nz_pcond

  !=============================================================================

  subroutine finalize_nz_pcond

    deallocate(izzftcd)
    deallocate(izcophc)
    deallocate(izcophg)
    deallocate(iztag1d)
    deallocate(ztpar  )

  end subroutine finalize_nz_pcond

  !=============================================================================

end module cs_nz_condensation
