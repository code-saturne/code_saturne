!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file cs_c_bindings.f90
!> Definition of structure, subroutine, and function C bindings.

module cs_c_bindings

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function checking the presence of a control file
    ! and dealing with the interactive control.

    subroutine cs_control_check_file()  &
      bind(C, name='cs_control_check_file')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_control_check_file

    !---------------------------------------------------------------------------

    !> \endcond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================


  !=============================================================================

end module cs_c_bindings
