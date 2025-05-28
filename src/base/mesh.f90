!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2025 EDF S.A.
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

!> \file mesh.f90
!> \brief Module for mesh-related arrays

module mesh

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup mesh Mesh Fortran structure, pointers to the C structure

  !> \addtogroup mesh
  !> \{

  !> \anchor ndim
  !> spatial dimension (3)
  integer :: ndim
  parameter(ndim=3)

  !> \anchor ncelet
  !> number of extended (real + ghost of the 'halo') cells. See \ref note_1
  integer, save :: ncelet = 0

  !> \anchor ncel
  !> number of real cells in the mesh
  integer, save :: ncel = 0

  !> \anchor nfabor
  !> number of boundary faces (see \ref note_2)
  integer, save :: nfabor = 0

  ! pointer to C array used by ifabor (0 to n-1 numbering)
  integer, dimension(:), pointer :: ifabor_0

  !> \anchor xyzcen
  !> coordinate of the cell centers
  double precision, dimension(:,:), pointer :: xyzcen

  !> \anchor cdgfbo
  !> coordinates of the centers of the boundary faces
  double precision, dimension(:,:), pointer :: cdgfbo

  !> \anchor surfbn
  !> norm of the surface of the boundary faces
  double precision, dimension(:), pointer :: surfbn

  !=============================================================================

contains

  !=============================================================================

  !> \anchor ifabor
  !> index-number of the (unique) neighboring cell for each boundary face

  elemental pure function ifabor(ifac) result(icel)

    implicit none

    ! Parameters

    integer, intent(in) :: ifac
    integer             :: icel

    ! Function body

    icel = ifabor_0(ifac) + 1

  end function ifabor

  !> \}

end module mesh
