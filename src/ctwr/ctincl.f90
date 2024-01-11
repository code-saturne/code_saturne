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

! Module for cooling towers

module ctincl

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use cstphy

  implicit none

  !=============================================================================

  !> \defgroup ctincl Module for cooling towers constants

  !> \addtogroup ctincl
  !> \{

  !=============================================================================

  !> Cp of dry air
  real(c_double), pointer, save :: cp_a

  !> Cp of water vapor
  real(c_double), pointer, save :: cp_v

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    ! Interface to C function retrieving pointers to members of the
    ! global fluid properties structure
    subroutine cs_air_glob_properties_get_pointers(cp_a, cp_v) &
        bind(C, name='cs_air_glob_properties_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: cp_a, cp_v
    end subroutine cs_air_glob_properties_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

contains

  subroutine ctwr_properties_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_cp_a, c_cp_v

    call cs_air_glob_properties_get_pointers(c_cp_a, c_cp_v)

    call c_f_pointer(c_cp_a, cp_a)
    call c_f_pointer(c_cp_v, cp_v)

  end subroutine ctwr_properties_init

  !=============================================================================

end module ctincl
