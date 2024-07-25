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

  !> Activates (> 0) or not (= 0), activate the ALE module:
  !>  - 0: not activated
  !>  - 1: legacy solver
  !>  - 2: CDO solver
  integer(c_int), pointer, save :: iale

  !> \}

contains

  !=============================================================================

  subroutine map_ale

    use, intrinsic :: iso_c_binding

    !---------------------------------------------------------------------------

    interface

      subroutine cs_f_ale_get_pointers(iale) &
        bind(C, name='cs_f_ale_get_pointers')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), intent(out) :: iale
      end subroutine cs_f_ale_get_pointers

    end interface

    ! Local variables

    type(c_ptr) :: c_iale

    call cs_f_ale_get_pointers(c_iale)

    call c_f_pointer(c_iale, iale)

  end subroutine map_ale

  !=============================================================================

end module albase
