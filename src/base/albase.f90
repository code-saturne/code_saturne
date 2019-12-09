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

!> \file albase.f90
!> Module for Arbitrary Lagrangian Eulerian method (ALE)

module albase

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \defgroup albase Module for Arbitrary Lagrangian Eulerian method (ALE)

  !> \addtogroup albase
  !> \{

  !> Activates (=1 o) or not (=0), activate the ALE module:
  !>  - 0: not activated
  !>  - 1: legacy solver
  !>  - 2: CDO solver
  integer(c_int), pointer, save :: iale
  !> the number of sub-iterations of initialization of the fluid
  integer, save :: nalinf
  !> maximum number of implicitation iterations of the structure displacement
  integer, save :: nalimx
  !> relative precision of implicitation of the structure displacement
  double precision, save :: epalim
  !> iteration (yes=1, no=0) to initialize ALE
  integer, save :: italin

  !> indicator of imposed displacement
  integer, allocatable, dimension(:) :: impale
  !> defines the mesh velocity from the color of the boundary faces,
  !> or more generally from their properties (colors, groups, ...),
  !> from the boundary conditions defined in cs user boundary conditions,
  !> or even from their coordinates.
  integer, allocatable, dimension(:) :: ialtyb
  !> initial mesh coordinates
  double precision, allocatable, dimension(:,:) :: xyzno0
  !> Pointer to field over vertices: mesh displacement
  integer, save :: fdiale
  !> \}

contains

  !=============================================================================

  subroutine init_ale (nfabor, nnod)

    use cplsat
    use optcal

    ! Arguments

    integer, intent(in) :: nfabor, nnod

    if (iale.ge.1) then
      allocate(xyzno0(3,nnod))
      allocate(impale(nnod))
      allocate(ialtyb(nfabor))
    endif

  end subroutine init_ale

  !=============================================================================

  subroutine finalize_ale

    use cplsat

    if (iale.ge.1) then
      deallocate(xyzno0)
      deallocate(impale)
      deallocate(ialtyb)
    endif

  end subroutine finalize_ale

  !=============================================================================

  subroutine map_ale

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    ! Local variables

    type(c_ptr) :: c_iale

    call cs_f_ale_get_pointers(c_iale)

    call c_f_pointer(c_iale, iale)

  end subroutine map_ale

  !=============================================================================

end module albase
