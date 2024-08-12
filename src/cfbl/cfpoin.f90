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

!> \file cfpoin.f90
!> Module for fuel combustion

module cfpoin

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \addtogroup compressible
  !> \{

  !> thermodynamic variables indicator for initialization
  !> mapping cs_cf_model_t::ithvar
  integer(c_int), pointer, save :: ithvar

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global compressible model structure

    subroutine cs_f_cf_model_get_pointers(ithvar)  &
      bind(C, name='cs_f_cf_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ithvar
    end subroutine cs_f_cf_model_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran compressible model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine cf_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ithvar

    call cs_f_cf_model_get_pointers(c_ithvar)

    call c_f_pointer(c_ithvar, ithvar)

  end subroutine cf_model_init

  !=============================================================================

end module cfpoin
