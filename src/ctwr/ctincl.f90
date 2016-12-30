!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

! Module for cooling towers

module ctincl

  !=============================================================================

  implicit none

  !=============================================================================


  ! cp_a     : Cp of dry air
  ! cp_v     : Cp of water vapour
  ! cp_l     : Cp of liquid water
  ! hv0      : enthalpy of vapourisation of water
  ! rho_l    : density of liquid water
  ! lambda_l : conductivity of liquid water
  ! lambda_h : conductivity of humid air

  double precision, save :: cp_a, cp_v, cp_l, hv0, rho_l, lambda_l, lambda_h

  ! humidity0 : initial absolute humidity in the cooling tower

  double precision, save :: humidity0

  !=============================================================================

end module ctincl
