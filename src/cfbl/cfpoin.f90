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

!> \file cfpoin.f90
!> Module for fuel combustion

module cfpoin

  !=============================================================================

  implicit none

  !=============================================================================

  ! ithvar : initialized thermodynamic variables indicator
  integer, save :: ithvar

  !=============================================================================

  ! Tableau Dimension       Description
  ! ifbet  ! nfabor        ! imposed thermal flux indicator at the boundary
  !                          (some boundary contributions of the total energy
  !                           equation have to be cancelled)
  ! icvfli ! nfabor        ! specific boundary convection flux indicator
  !                          for a Rusanov or an analytical flux
  !                          (some boundary contributions of the momentum
  !                           equation have to be cancelled)

  integer, allocatable, dimension(:) :: ifbet , icvfli

contains

  !=============================================================================

  subroutine init_compf ( nfabor)

    implicit none

    integer nfabor

    allocate(ifbet(nfabor))
    allocate(icvfli(nfabor))

  end subroutine init_compf

  !=============================================================================

  subroutine finalize_compf

    use ppincl
    implicit none

    deallocate(ifbet, icvfli)

  end subroutine finalize_compf

end module cfpoin
