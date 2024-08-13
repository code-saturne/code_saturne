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

!> \file cdomod.f90
!> \brief Store the mode of activation of CDO-HHO schemes
!>
!------------------------------------------------------------------------------

module cdomod

  !===========================================================================

  use, intrinsic :: iso_c_binding

  !=============================================================================

  implicit none

  !=============================================================================

  !> Activated (=1 or =2) or not activated (=0)
  !> If icdo=1 (CDO and FV at the same time)
  !> If icdo=2 (CDO only)
  integer(c_int), pointer, save :: icdo

  interface

    ! Interface to C function retrieving pointers to global CDO parameters

    subroutine cs_f_cdo_get_pointers(icdo) &
      bind(C, name='cs_f_cdo_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: icdo
    end subroutine cs_f_cdo_get_pointers

 end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran CDO flag.
  !> This maps Fortran pointers to global C structure members.

  subroutine cdo_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_icdo

    call cs_f_cdo_get_pointers(c_icdo)

    call c_f_pointer(c_icdo, icdo)

  end subroutine cdo_init

  !=============================================================================

end module cdomod
