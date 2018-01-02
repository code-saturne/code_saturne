!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file period.f90
!> \brief Module for periodicity flags

module period

  !=============================================================================

  implicit none

  !> \defgroup period Module for periodicity flags

  !> \addtogroup period
  !> \{

  !> presence of periodicity.
  !> - 1: periodicity is activated
  !> - 0: no periodicity (default value)
  integer, save :: iperio

  !> number of rotation periodicities. automaticly evaluated.
  !> default value is 0
  integer, save :: iperot

  ! TODO
  integer, save :: igrper

  !=============================================================================

  !> \}

end module period
