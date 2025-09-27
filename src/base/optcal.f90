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

!> \file optcal.f90
!> Module for calculation options

module optcal

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \defgroup optcal Module for calculation options

  !> \addtogroup optcal
  !> \{

  !----------------------------------------------------------------------------
  ! Space discretisation
  !----------------------------------------------------------------------------

  !> Indicator of a calculation restart (=1) or not (=0).
  !> This value is set automatically by the code; depending on
  !> whether a restart directory is present, and should not be modified by
  !> the user (no need for C mapping).
  integer, save :: isuite = 0

  !----------------------------------------------------------------------------
  ! Time stepping options
  !----------------------------------------------------------------------------

  !> \defgroup time_step_options Time step options and variables

  !> \addtogroup time_step_options
  !> \{

  !> Current absolute time step number.
  !> In case of restart, this is equal to ntpabs + number of new iterations.
  integer(c_int), pointer, save :: ntcabs

  !> Absolute time value for previous calculation.
  !>
  !> Current absolute time.
  !>
  !> For the restart calculations, \ref ttcabs takes
  !> into account the physical time of the previous calculations.\n
  !> If the time step is uniform (\ref idtvar = 0 or 1), \ref ttcabs
  !> increases of \ref dt (value of the time step) at each iteration.
  !> If the time step is non-uniform (\ref idtvar=2), \ref ttcabs
  !> increases of \ref dtref at each time step.\n
  !> \ref ttcabs} is initialised and updated automatically by the code,
  !> its value is not to be modified by the user.
  real(c_double), pointer, save :: ttcabs

  !> option for a variable time step
  !>    - -1: steady algorithm
  !>    -  0: constant time step
  !>    -  1: time step constant in space but variable in time
  !>    -  2: variable time step in space and in time
  !> If the numerical scheme is a second-order in time, only the
  !> option 0 is allowed.
  integer(c_int), pointer, save :: idtvar

  !> Reference time step
  !>
  !> This is the time step value used in the case of a calculation run with a
  !> uniform and constant time step, i.e. \ref idtvar =0 (restart calculation
  !> or not). It is the value used to initialize the time step in the case of
  !> an initial calculation run with a non-constant time step(\ref idtvar=1 or
  !> 2). It is also the value used to initialise the time step in the case of
  !> a restart calculation in which the type of time step has been changed
  !> (for instance, \ref idtvar=1 in the new calculation and \ref idtvar = 0 or
  !> 2 in the previous calculation).\n
  !> See \subpage user_initialization_time_step for examples.
  real(c_double), pointer, save :: dtref

  !> \}

  !----------------------------------------------------------------------------
  ! Transported scalars parameters
  !----------------------------------------------------------------------------

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global time step structure

    subroutine cs_f_time_step_get_pointers(nt_cur, dt_ref, t_cur)   &
      bind(C, name='cs_f_time_step_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nt_cur
      type(c_ptr), intent(out) :: dt_ref, t_cur
    end subroutine cs_f_time_step_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global time step options structure

    subroutine cs_f_time_step_options_get_pointers(idtvar)   &
      bind(C, name='cs_f_time_step_options_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: idtvar
    end subroutine cs_f_time_step_options_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize isuite

  subroutine indsui () &
    bind(C, name='cs_f_indsui')

    use, intrinsic :: iso_c_binding
    use cs_c_bindings
    implicit none

    isuite = cs_restart_present()

  end subroutine indsui

  !> \brief Initialize Fortran time step API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ntcabs
    type(c_ptr) :: c_dtref, c_ttcabs

    call cs_f_time_step_get_pointers(c_ntcabs, c_dtref, c_ttcabs)

    call c_f_pointer(c_ntcabs, ntcabs)
    call c_f_pointer(c_dtref,  dtref)
    call c_f_pointer(c_ttcabs, ttcabs)

  end subroutine time_step_init

  !> \brief Initialize Fortran time step options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_idtvar

    call cs_f_time_step_options_get_pointers(c_idtvar)

    call c_f_pointer(c_idtvar, idtvar)

  end subroutine time_step_options_init

  !=============================================================================

end module optcal
