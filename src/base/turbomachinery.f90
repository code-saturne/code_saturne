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

!> \file turbomachinery.f90
!> \brief Module for turbomachinery computations

module turbomachinery

  !=============================================================================
  !> \defgroup at_turbomachinery Module for turbomachinery computations

  !> \addtogroup at_turbomachinery
  !> \{

  !> Type of turbomachinery computation:
  !>   none (0), frozen rotor (1), transient (2)

  integer, save :: iturbo

  !> \}

  !=============================================================================

  interface

    ! Interface to C function mapping some data for turbomachinery

    subroutine map_turbomachinery_model(iturbo2) &
      bind(C, name='cs_f_map_turbomachinery_model')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(out) :: iturbo2
    end subroutine map_turbomachinery_model

  end interface

  !=============================================================================

end module turbomachinery
