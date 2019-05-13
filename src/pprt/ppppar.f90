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

!> \file ppppar.f90
!> General module for specific physics containing common parameters

module ppppar

  !=============================================================================

  implicit none

  !===========================================================================


  !> \defgroup ppppar General module for specific physics containing common parameters

  !> \addtogroup ppppar
  !> \{

  !> maximum number of boundary zones
  integer    nbzppm
  parameter (nbzppm=2000)

  !> maximum index of boundary zones
  integer    nozppm
  parameter (nozppm=2000)

  !--> Pointeurs variables combustion charbon pulverise cpincl, ppincl

  !    ncharm --> nombre maximal de charbons
  !    ncpcmx --> nombre maximal de classes par charbon
  !    nclcpm --> Nombre total de classes

  !> maximum number of coals
  integer    ncharm
  !> maximum number of coals classes
  integer    ncpcmx

  !> maximum number of coal classes for the pulverised coal combustion module
  integer    nclcpm
  parameter (ncharm=5, ncpcmx=20, nclcpm=ncharm*ncpcmx)

  !=============================================================================

  !> \}

end module ppppar
