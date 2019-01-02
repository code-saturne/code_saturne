!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_modules-user-arrays.f90
!>
!> \brief User-defined module example
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

module user_module

  !=============================================================================

  implicit none

  !=============================================================================
!< [variables]
  ! Example: allocatable user arrays

  integer,          dimension(:), allocatable :: iwork
  double precision, dimension(:,:), pointer :: rwork => null()
!< [variables]
contains

  !=============================================================================
!< [allocate]
  ! Allocate arrays

  subroutine init_user_module(ncel, ncelet)

    ! Arguments

    integer, intent(in) :: ncel, ncelet

    ! Local variables

    integer :: err = 0

    if (.not.allocated(iwork)) then
      allocate(iwork(ncelet), stat=err)
    endif

    if (err .eq. 0 .and. .not.associated(rwork)) then
      allocate(rwork(3, ncelet), stat=err)
    endif

    if (err /= 0) then
      write (*, *) "Error allocating array."
      call csexit(err)
    endif

    return

  end subroutine init_user_module
!< [allocate]
  !=============================================================================
!< [c_pointer]
  ! Pass pointer to rwork array to C

  function get_user_module_rwork() result(r) &
    bind(C, name='get_user_module_rwork')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr) :: r
    if (associated(rwork)) then
      r = c_loc(rwork(1,1))
    else
      r = c_null_ptr
    endif
  end function get_user_module_rwork
!< [c_pointer]
  !=============================================================================
!< [free]
  ! Free related arrays

  subroutine finalize_user_module

    if (allocated(iwork)) then
      deallocate(iwork)
    endif

    if (associated(rwork)) then
      deallocate(rwork)
    endif

  end subroutine finalize_user_module
!< [free]

  !=============================================================================

end module user_module

!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
