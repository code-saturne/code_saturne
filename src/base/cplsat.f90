!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cplsat.f90
!> \brief Module for code/code coupling

module cplsat

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup cplsat Module for code/code coupling
  !> code / code - management of key parameters

  !> \addtogroup cplsat
  !> \{

  !> number of couplings code_saturne / code_saturne
  integer, save :: nbrcpl
  !> indicator coupling face / face only
  integer, save :: ifaccp
  !> maximum permissible number of coupling
  integer   nbcpmx
  parameter(nbcpmx=10)
  !> turbulence model of the remote instance
  integer, save :: iturcp(nbcpmx)
  !> indicator to update location of the coupling
  integer, save :: imajcp(nbcpmx)
  !> indicator of calulation in relative reference frame
  integer, save :: icormx(nbcpmx)
  !> number of variables to send/receive
  integer, save :: nvarcp(nbcpmx)
  !> size of exchange tables
  integer, save :: nvarto(nbcpmx)
  !> Absolute time value after the mesh starts to rotate (if it does),
  !> for previous calculation
  double precision, save :: ttpmob
  !> Current absolute time after the mesh starts to rotate (if it does).
  !> In case of restart, this is equal to ttpmob + additional computed time.
  double precision, save :: ttcmob
  !> \}

  !=============================================================================

end module cplsat
