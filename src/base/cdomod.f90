!-------------------------------------------------------------------------------

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

  !> Activated (=1 or =2) or not (=0)
  !> If icdo=1 (CDO and FV)
  !> If icdo=2 (CDO only)
  integer, save :: icdo

end module cdomod

!=============================================================================

subroutine set_cdo_mode (icdoval)

  !===========================================================================
  ! Module files
  !===========================================================================

  use cdomod

  !===========================================================================

  implicit none

  ! Arguments

  integer(c_int)          icdoval

  !===========================================================================

  icdo = icdoval

  return
end subroutine set_cdo_mode
