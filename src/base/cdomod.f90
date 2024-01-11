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

    ! Interface to C function to solve the unsteady state for related CDO
    ! equations

    subroutine cs_f_cdo_solve_unsteady_state_domain()  &
      bind(C, name='cs_f_cdo_solve_unsteady_state_domain')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_cdo_solve_unsteady_state_domain

    ! Interface to C function to solve the steady state for related CDO
    ! equations

    subroutine cs_f_cdo_solve_steady_state_domain()  &
      bind(C, name='cs_f_cdo_solve_steady_state_domain')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_cdo_solve_steady_state_domain

    ! Interface to C function related to the initialization CDO systems

    subroutine cs_f_domain_initialize_cdo_systems()  &
      bind(C, name='cs_f_domain_initialize_cdo_systems')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_domain_initialize_cdo_systems

    ! Interface to C function to postprocess data related to CDO schemes

    subroutine cs_f_cdo_post_domain()  &
      bind(C, name='cs_f_cdo_post_domain')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_cdo_post_domain

    ! Interface to C function to force the resolution of steady equation

    subroutine cs_equation_solve_steady_state_wrapper(eqname) &
      bind(C, name='cs_equation_solve_steady_state_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: eqname
    end subroutine cs_equation_solve_steady_state_wrapper

    ! Interface to C function to force the resolution of an unsteady equation

    subroutine solve_cdo_equation(cur2prev, eqname) &
      bind(C, name='cs_equation_solve_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(c_bool), value :: cur2prev
      character(kind=c_char, len=1), dimension(*), intent(in) :: eqname
    end subroutine solve_cdo_equation

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

  subroutine solve_steady_state_cdo_equation(eqname)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: eqname

    ! Local variables

    character(len=len_trim(eqname)+1, kind=c_char) :: c_eqname

    c_eqname = trim(eqname)//c_null_char

    call cs_equation_solve_steady_state_wrapper(c_eqname)

    return

  end subroutine solve_steady_state_cdo_equation

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

