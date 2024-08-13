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

!> \file vof.f90
!> \brief Module for Volume-Of-Fluid method

module vof

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  interface

     !---------------------------------------------------------------------------

     !> \cond DOXYGEN_SHOULD_SKIP_THIS

     !---------------------------------------------------------------------------

     ! Interface to C function retrieving pointers to VOF model indicator
     ! and parameters

     subroutine cs_f_vof_get_pointers(ivofmt) &
       bind(C, name='cs_f_vof_get_pointers')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(out) :: ivofmt
     end subroutine cs_f_vof_get_pointers

     !---------------------------------------------------------------------------

     !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

     !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran VOF model API.
  !> This maps Fortran pointers to global C structure members and indicator.

  subroutine vof_model_init

    use, intrinsic :: iso_c_binding
    use optcal, only:ivofmt

    implicit none

    ! Local variables

    type(c_ptr) :: c_ivofmt

    call cs_f_vof_get_pointers(c_ivofmt)

    call c_f_pointer(c_ivofmt, ivofmt)

  end subroutine vof_model_init

end module vof
